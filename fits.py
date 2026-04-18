"""
@author: micha
"""

import numpy as np
from scipy.optimize import curve_fit
import convert_data
import fit_models 

# substract regions
def exclude_region(x, y, x_min, x_max):
    idx = np.where(np.logical_or(x<=x_min, x>=x_max))
    return x[idx], y[idx]

def DD_line(data, exclude=[[None, None]], model_types=['linear'], inits=[dict({'k':1,'d':1e-8})], sigma_lor=200, sigma_baskov=200, bounds=[None], **kwargs):
    # prepare data
    x = data[:,0]
    y = data[:,1]
    if exclude[0][0] is not None:
        for i in range(len(exclude)):
            x, y = exclude_region(x,y,exclude[i][0],exclude[i][1])     
    if len(bounds) < len(model_types):
        bounds.extend([None for i in range(len(model_types)-len(bounds))])
    y_stretch = np.max(y)
    y = convert_data.max_to_one(y)
    # define fit model
    comp_model = fit_models.Model()
    for model, init in zip(model_types, inits):
        if model == 'voigt':
            comp_model.add_fit_func(model,init,bounds=fit_models.bounds_voigt())
        elif model == 'lorentzian':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-40,0], [np.inf,init['x0']+40,sigma_lor]])#80
        elif model == 'baskovian':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-40,0], [np.inf,init['x0']+40,sigma_baskov]])# 40 = Presentation
        else:
            comp_model.add_fit_func(model,init)
    res, cov = curve_fit(comp_model.get_fit_func(), x, y, p0=comp_model.get_params(), bounds=comp_model.get_bounds())
    print(res,y_stretch)
    #comp_model.set_fit_params(res, y_stretch=y_stretch)
    comp_model.set_fit_params(res, y_stretch=1)
    return comp_model, res

def G_line(data, exclude=[[None, None]], model_types=['linear'], inits=[dict({'k':1,'d':1e-8})], bounds=[None], **kwargs):
    # prepare data
    x = data[:,0]
    y = data[:,1]
    if exclude[0][0] is not None:
        for i in range(len(exclude)):
            x, y = exclude_region(x,y,exclude[i][0],exclude[i][1])     
    if len(bounds) < len(model_types):
        bounds.extend([None for i in range(len(model_types)-len(bounds))])
    if 'normalize' in kwargs and kwargs['normalize'] is True:
        y_stretch = 1
    else:
        y_stretch = np.max(y)
    y = convert_data.max_to_one(y)
    # define fit model
    comp_model = fit_models.Model()
    for model, init in zip(model_types, inits):
        if model == 'lorentzian':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-10,0], [np.inf,init['x0']+10,np.inf]])#fit_models.bounds_lorentzian())
        elif model == 'voigt':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-10,0,0], [np.inf,init['x0']+10,10,np.inf]])
        elif model == 'bwf':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-10,0,-np.inf], [np.inf,init['x0']+10,np.inf,0]])
        elif model == 'polynomial':
            comp_model.add_fit_func(model,init,bounds=[[0],[np.inf]])
        else:
            comp_model.add_fit_func(model,init)
        #if model == 'voigt':
        #    comp_model.add_fit_func(model,init,bounds=fit_models.bounds_voigt())
        #else:
        #    comp_model.add_fit_func(model,init)
    comp_model.renorm_fits()
    res, cov = curve_fit(comp_model.get_fit_func(), x, y, p0=comp_model.get_params(), bounds=comp_model.get_bounds())
    comp_model.set_fit_params(res, y_stretch=y_stretch)
    
    
    def GFG(arr):
        np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
        print(arr)
    temp_print = comp_model.get_params(separate=True)
    for p in temp_print:
        GFG(np.array(p))
    print(len(temp_print))
    
    return comp_model, res

def RBM_lines(data, exclude=[[None, None]], model_types=['linear'], inits=[dict({'k':1,'d':1e-8})], sigma=50, **kwargs):
    # prepare data
    x = data[:,0]
    y = data[:,1]
    if exclude[0][0] is not None:
        for i in range(len(exclude)):
            x, y = exclude_region(x,y,exclude[i][0],exclude[i][1])   
    x_min = 80
    x_max = 300
        
    y = convert_data.max_to_one(y)
    if 'normalize' in kwargs and kwargs['normalize'] is True:
        y_stretch = 1
    else:
        y_stretch = np.max(y)
    y = convert_data.max_to_one(y)
    # define fit model
    comp_model = fit_models.Model()
    for model, init in zip(model_types, inits):
        if model == 'voigt':
            comp_model.add_fit_func(model,init,bounds=fit_models.bounds_voigt())
        elif model == 'polynomial':
            comp_model.add_fit_func(model,init,bounds=[[0],[np.inf]])
        elif model == 'lorentzian':
            comp_model.add_fit_func(model,init,bounds=[[0,init['x0']-10,1], [np.inf,init['x0']+10,sigma]])#fit_models.bounds_lorentzian())
        elif model == 'baskovian':
            comp_model.add_fit_func(model,init,bounds=fit_models.bounds_baskovian())
        elif model == 'gaussian':
            comp_model.add_fit_func(model,init,bounds=[[0,0,0],[np.inf,100,np.inf]])
        else:
            comp_model.add_fit_func(model,init)
    comp_model.renorm_fits()
    print(comp_model.get_params())
    res, cov = curve_fit(comp_model.get_fit_func(), x, y, p0=comp_model.get_params(), bounds=comp_model.get_bounds())
    comp_model.set_fit_params(res, y_stretch=y_stretch)
    
    def GFG(arr):
        np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
        print(arr)
    temp_print = comp_model.get_params(separate=True)
    for p in temp_print:
        GFG(np.array(p))
    print(len(temp_print))
    #print(comp_model.get_params(separate=True))
    
    return comp_model, res, x,y
