"""
@author: Michael Labenbacher
"""
import numpy as np
from scipy.optimize import curve_fit
import convert_data
import fit_models

def normalize_counts(y):
    return y/np.max(y)

def gauss(x, c, a, x0, sigma):
    return c + a*np.exp(-(x-x0)**2/(2*sigma**2))# XXXXX: Use fit_models gauss function instead !!!

# substract regions
def exclude_region(x, y, x_min, x_max):
    idx = np.where(np.logical_or(x<=x_min, x>=x_max))
    return x[idx], y[idx]

def correction_center_rayleigh(data, label='633', xmin=-10., xmax=10., p0=[0.1, 1., 0., 5.], **kwargs):
    x = np.array([d[0] for d in data], dtype=np.float64)
    y = np.array([d[1] for d in data], dtype=np.float64)
    y_fit = normalize_counts(y)
    idx = np.where(np.logical_and(x<xmax, x>xmin))
    x_fit = x[idx]
    y_fit = y_fit[idx]
    # XXXXX: Define a model and use fit guass there with 0th order polynom
    res, cov = curve_fit(gauss, x_fit, y_fit, p0=p0)
    if label == '405':
        x_res = x-(res[2]-kwargs['wavelength'])
    else:
        x_res = x-res[2]
    y_res = y-res[0]
    return x_res, y_res

def correction_background(data, exclude=[[None, None]], model_types=['linear'], init=[dict({'k':1,'d':1e-8})], bounds=[None], **kwargs):
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
    print(comp_model.get_fit_params())
    res, cov = curve_fit(comp_model.get_fit_func(), x, y, p0=comp_model.get_fit_params(), bounds=comp_model.get_bounds())
    #res, cov = curve_fit(func_poly, x, y, p0=comp_model.get_fit_params())
    print(res)
    comp_model.set_fit_params(res, y_stretch=temp_ymax)

    if 'res' in kwargs:
        return x, y*temp_ymax, comp_model, res, cov
    else: 
        return x, y*temp_ymax, comp_model
    
    