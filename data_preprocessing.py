"""
@author: Michael Labenbacher
"""
import numpy as np
from scipy.optimize import curve_fit
import convert_data
import fit_models
import plots

def get_wavelength_from_label(data_key, file_number=True):
    key_sep = data_key.split('-')
    if file_number is True:
        if len(key_sep)>2:
            return key_sep[0], int(key_sep[2])
        else: 
            return key_sep[0], int(key_sep[1])
    else:
        return key_sep[0]

def get_filequantity_for_wavelength(data_keys, wavelength='633'):
    quantity = 0 
    for key in data_keys:
        key_wavelength = key.split('-')[0]
        if key_wavelength == wavelength:
            quantity += 1
    return quantity

def get_labels_for_wavelength(data_keys, wavelength='633'):
    labels = []
    for key in data_keys:
        key_sep = key.split('-')
        if key_sep[0] == wavelength:
            labels.append(key)
    return labels

def get_datas_for_wavelength(data, wavelength='633'):
    labels = get_labels_for_wavelength(data.keys(), wavelength)
    data_subset = {l: data[l] for l in labels}
    return data_subset

def normalize_counts(y):
    return y/np.max(y)

# substract regions
def exclude_region(x, y, x_min, x_max):
    idx = np.where(np.logical_or(x<=x_min, x>=x_max))
    return x[idx], y[idx]

def correction_center_rayleigh(data, label='633', xmin=-10., xmax=10., p0=[0.1, 1., 0., 5.], plot=False, **kwargs):
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    y_fit = normalize_counts(y)
    idx = np.where(np.logical_and(x<xmax, x>xmin))
    x_fit = x[idx]
    y_fit = y_fit[idx]
    
    comp_model = fit_models.Model()
    comp_model.add_fit_func('polynomial',{'a0':p0[0]})
    comp_model.add_fit_func('gaussian',{'c':p0[1],'x0':p0[2],'sigma':p0[3]})
    res, cov = curve_fit(comp_model.get_fit_func(), x_fit, y_fit, p0=comp_model.get_params(), bounds=comp_model.get_bounds())
    comp_model.set_fit_params(res)

    x_res_fit = np.linspace(xmin, xmax, num=1000, endpoint=True)
    if label in get_labels_for_wavelength([label],'405'):
        x_res = x-(res[2]-kwargs['wavelength'])# offset x from gaussian x0 in nm (correct difference to 405nm)
        if plot is True:
            x_res_fit += -(res[2]-kwargs['wavelength'])
    else:
        x_res = x-res[2]# offset x from gaussian x0 in cm-1
        if plot is True:
            x_res_fit += -res[2]
    y_res = y-res[0]# offset y from polynomial a0
    if plot is True:
        x_fit_ = np.linspace(xmin, xmax, num=1000, endpoint=True)
        plots.rayleigh(init=(x_fit,y_fit), 
                       init_fit=(x_fit_,comp_model.get_fit_func()(x_fit_,*comp_model.get_params())), 
                       res=(x_res[idx],y_res[idx]),
                       res_fit=(x_res_fit,comp_model.get_fit_func()(x_fit_,*comp_model.get_params())), 
                       plot_range=(xmin,xmax), label=label)
        print(res)
    return x_res, y_res

def correction_background(data, exclude=[[None, None]], epsilon=1e-6):
    # prepare data
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    temp_ymax = np.max(y)
    y = convert_data.max_to_one(y)
    x_ = np.copy(x)
    y_ = np.copy(y)
    if exclude[0][0] is not None:
        for i in range(len(exclude)):
            x_, y_ = exclude_region(x_,y_,exclude[i][0],exclude[i][1])     
    
    y = y - np.min(y_) + epsilon
    return x,y*temp_ymax

def correction_background_(data, exclude=[[None, None]], model_types=['linear'], init=[dict({'k':1,'d':1e-8})], bounds=[None], **kwargs):
    # prepare data
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    temp_ymax = np.max(y)
    y = convert_data.max_to_one(y)
    if exclude[0][0] is not None:
        for i in range(len(exclude)):
            x, y = exclude_region(x,y,exclude[i][0],exclude[i][1])     
    if len(bounds) < len(model_types):
        bounds.extend([None for i in range(len(model_types)-len(bounds))])
    # generate model
    comp_model = fit_models.Model()
    for i, model_type in enumerate(model_types):   
        if model_type in ['linear']:
            #comp_model.params_update(list(init[i].values()))
            if all(key in init[i] for key in ('k', 'd')):
                init_list = [init[i][key] for key in ['k', 'd']]
                comp_model.params_update(init_list)
            else:
                print("Linear model needs initial values k and d.")
            comp_model.funcs_update(model_type)
        elif model_type in ['polynomial']:
            init_list = []
            for a_i in range(len(init[i])):
                if 'a'+str(a_i) in init[i]:
                    init_list.append(init[i]['a'+str(a_i)])
                else:
                    print("Polynomial model only supports a_i with i=0,1,2... as initial values.")
                    break 
            comp_model.params_update(init_list)
            comp_model.funcs_update(model_type)
        elif model_type in ['gaussian']:
            if all(key in init[i] for key in ('c', 'x0', 'sigma')):
                init_list = [init[i][key] for key in ['c', 'x0', 'sigma']]
                comp_model.params_update(init_list, bounds=fit_models.bounds_gaussian())
            else:
                print("Gaussian model needs initial values c, x0 and sigma.")
            comp_model.funcs_update(model_type)
        elif model_type in ['lorentzian']:
            keys_required = [key in init[i] for key in ('c', 'x0', 'gamma')]
            if all(keys_required):
                init_list = [None]*3
                for i_k,key in enumerate(['c', 'x0', 'gamma']):
                    init_list[i_k] = init[i][key]
                    if key == 'c':
                        init_list[i_k] *= 1/fit_models.func_lorentzian(init[i]['x0'], 1, init[i]['x0'], init[i]['gamma'])
                if bounds[i] is not None:
                    bounds_add = bounds[i]
                else: 
                    bounds_add = fit_models.bounds_lorentzian()
                comp_model.params_update(init_list, bounds=bounds_add)
            else:
                print("Lorentzian model needs initial values c, x0 and gamma.")
            if len(init[i].keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
            comp_model.funcs_update(model_type)
        elif model_type in ['voigt']:
            if all(key in init[i] for key in ('c', 'x0', 'sigma', 'gamma')):
                init_list = [None]*4
                for i_k,key in enumerate(['c', 'x0', 'sigma', 'gamma']):
                    init_list[i_k] = init[i][key]
                    if key == 'c':
                        init_list[i_k] *= 1/fit_models.func_voigt(init[i]['x0'], 1, init[i]['x0'], init[i]['sigma'], init[i]['gamma'])
                #init_list = [init[i][key] for key in ['c', 'x0', 'sigma', 'gamma']]
                if bounds[i] is not None:
                    bounds_add = bounds[i]
                else: 
                    bounds_add = fit_models.bounds_voigt()
                comp_model.params_update(init_list, bounds=bounds_add)
            else:
                print("Voigt model needs initial values c, x0, sigma and gamma.")
            comp_model.funcs_update(model_type)
        else:
            print("This model type {0} is not supported.".format(model_type))
    # fit data
    res, cov = curve_fit(comp_model.get_fit_func(), x, y, p0=comp_model.get_params(), bounds=comp_model.get_bounds())
    #res, cov = curve_fit(func_poly, x, y, p0=comp_model.get_fit_params())
    print(res)
    comp_model.set_fit_params(res, y_stretch=temp_ymax)

    if 'res' in kwargs:
        return x, y*temp_ymax, comp_model, res, cov
    else: 
        return x, y*temp_ymax, comp_model
    
def correction_spectrometer_y_init(data,pixels=20):
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    
    temp_a = np.mean(y[:pixels])
    y -= temp_a
    return np.column_stack((x, y))  
 
def correction_spectrometer_y_jump(data,center,pixels=10):
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    
    temp_a = np.mean(y[center:(center+pixels)])
    temp_b = np.mean(y[(center-pixels):center])
    y[center:] *= temp_b / temp_a
    
    return np.column_stack((x, y))   

def correction_spectrometer_outlier(data, smooth_points=[None], smooth_range=[None]):
    # prepare data
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    
    if smooth_points is not None:
        for s,sr in zip(smooth_points,smooth_range):
            y[s] = np.mean(y[s-sr:s+sr+1])
    
    return np.column_stack((x, y))   
    