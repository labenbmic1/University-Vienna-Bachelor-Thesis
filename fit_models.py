"""
@author: micha
"""

import numpy as np
from scipy.special import voigt_profile

def func_baskovian(x, *params):
    """ 3
    Baskovian fit function (3 parameters): 
    
    :param x: 1D array (x-data)
    :param params: c, x0, gamma (c=area, gamma=FWHM, gamma_=beta**2)
    
    :return: 1D array (y-data)
    """
    c, x0, gamma = params
    gamma_ = (gamma**2)/(4.*(2**(2/3)-1.))
    return c*gamma_**(3/2)/(2.* ( (x-x0)**2 + gamma_ )**(3/2) )
    #return c*gamma_/(2.* ( (x-x0)**2 + gamma_ )**(3/2) )
def func_baskovian_par(*params):
    """
    :return: x0, peak, FWHM
    """
    c, x0, gamma = params
    return func_baskovian(x0, *params), x0, gamma
def func_gaussian(x, *params):
    """ 3
    Gaussian fit function (3 parameters): c*exp(-(x-x0)**2/(2*sigma**2))
    
    :param x: 1D array (x-data)
    :param params: c, x0, sigma
    
    :return: 1D array (y-data)
    """
    return params[0]*np.exp(-(x-params[1])**2/(2*params[2]**2))
def func_lorentzian(x, *params):
    """ 3
    Lorentzian fit function (3 parameters): 
    
    :param x: 1D array (x-data)
    :param params: c, x0, gamma
    
    :return: 1D array (y-data)
    """
    c, x0, gamma = params
    return c * gamma**2 / ((x - x0)**2 + gamma**2)
def func_bwf(x, *params):
    """ 3
    Breit-Wigner Fanon fit function (3 parameters): 
    
    :param x: 1D array (x-data)
    :param params: c, x0, gamma, q
    
    :return: 1D array (y-data)
    """
    c, x0, gamma, q = params
    s = (x-x0)/gamma
    return c * (1+s/q)**2 / (1 + (s)**2)
def func_voigt(x, *params):
    """ 4
    Voigt fit function (4 parameters).
    
    :param x: 1D array (x-data)
    :param params: c, x0, sigma, gamma
    
    :return: 1D array (y-data)
    """
    return params[0]*voigt_profile(x-params[1], params[2], params[3])
def func_lin(x, *params):
    """ 2
    Linear fit function (2 parameters): k*x+d
    
    :param x: 1D array (x-data)
    :param params: k, d
    
    :return: 1D array (y-data)
    """
    return params[0]*x + params[1]
def func_poly(x, *params):
    """ n
    Polynomial fit function (n parameters): a0 + a1*x +...+ an*x**n
    
    :param x: 1D array (x-data)
    :param params: a0, a1,... , an
    
    :return: 1D array (y-data)
    """
    #return sum([p*(x**i) for i, p in enumerate(params)])
    return np.polynomial.Polynomial(params)(x)

# bounds or here also possible with abs in func definitions
def bounds_gaussian():
    return [0,-np.inf,0], [np.inf,np.inf,np.inf]
def bounds_voigt():
    return [0,-np.inf,0,0], [np.inf,np.inf,np.inf,np.inf]
def bounds_baskovian():
    return [0,0,0], [np.inf,np.inf,np.inf]
def bounds_lorentzian():
    return [0,0,0], [np.inf,np.inf,np.inf]
def bounds_bwf():
    return [0,0,0,-np.inf], [np.inf,np.inf,np.inf,0]


class Model(object):
    def __init__(self, model_type=None, params=[]):
        self._total = 0 
        self._num_params = []
        self._name_params = []
        self._params = []
        self._bounds_low = []
        self._bounds_high = []
        self._funcs = []
        self._model_types = []
        if model_type is not None:
            self._total += 1      
            self._num_params += [len(params)]
            self._params += [params] 
            self.funcs_update(model_type)
            if model_type in ['gaussian']:
                self._bounds_low, self._bounds_high = bounds_gaussian()
            elif model_type in ['voigt']:
                self._bounds_low, self._bounds_high = bounds_voigt()
            elif model_type in ['baskovian']:
                self._bounds_low, self._bounds_high = bounds_baskovian()
            elif model_type in ['bwf']:
                self._bounds_low, self._bounds_high = bounds_bwf()

    def params_update(self, params:list, bounds=None):
        self._params += [params]
        self._num_params += [len(params)]
        if bounds is None:
            bounds = [[-np.inf]*len(params), [np.inf]*len(params)]
        self.bounds_update(bounds[0], bounds[1])
    
    def bounds_update(self, bounds_low, bounds_high):
        self._bounds_low += [bounds_low]
        self._bounds_high += [bounds_high]
    
    def param_names_update(self, model_type, params:list):
        if model_type in ['linear']:
            self._name_params.extend(['k','d'])
        elif model_type in ['polynomial']:
            self._name_params.extend(['a'+str(i) for i in range(len(params))])
        elif model_type in ['gaussian']:
            self._name_params.extend(['c','x0','sigma'])
        elif model_type in ['lorentzian']:
            self._name_params.extend(['c','x0','gamma'])
        elif model_type in ['voigt']:
            self._name_params.extend(['c', 'x0', 'sigma', 'gamma'])
        elif model_type in ['baskovian']:
            self._name_params.extend(['c','x0','gamma'])
        elif model_type in ['bwf']:
            self._name_params.extend(['c','x0','gamma','q'])

    def funcs_update(self, model_type):
        if model_type in ['linear']:
            self._funcs += [func_lin]
            self._model_types += ['linear']
        elif model_type in ['polynomial']:
            self._funcs += [func_poly]
            self._model_types += ['polynomial']
        elif model_type in ['gaussian']:
            self._funcs += [func_gaussian]
            self._model_types += ['gaussian']
        elif model_type in ['lorentzian']:
            self._funcs += [func_lorentzian]
            self._model_types += ['lorentzian']
        elif model_type in ['voigt']:
            self._funcs += [func_voigt]
            self._model_types += ['voigt']
        elif model_type in ['baskovian']:
            self._funcs += [func_baskovian]
            self._model_types += ['baskovian']
        elif model_type in ['bwf']:
            self._funcs += [func_bwf]
            self._model_types += ['bwf']
    
    def set_fit_params(self, params, y_stretch=1):
        par_start = 0
        for i, (num_par, model_type) in enumerate(zip(self._num_params, self._model_types)):
            par_end = par_start + num_par
            if model_type in ['linear']:
                self._params[i] = [y_stretch*params[i] for i in range(par_start, par_end)]
            elif model_type in ['polynomial']:
                self._params[i] = [y_stretch*params[i] for i in range(par_start, par_end)]
            elif model_type in ['gaussian']:
                self._params[i] = [params[i] for i in range(par_start, par_end)]
                self._params[i][0] *= y_stretch
            elif model_type in ['lorentzian']:
                self._params[i] = [params[i] for i in range(par_start, par_end)]
                self._params[i][0] *= y_stretch
            elif model_type in ['voigt']:
                self._params[i] = [params[i] for i in range(par_start, par_end)]
                self._params[i][0] *= y_stretch
            elif model_type in ['baskovian']:
                self._params[i] = [params[i] for i in range(par_start, par_end)]
                self._params[i][0] *= y_stretch 
            elif model_type in ['bwf']:
                self._params[i] = [params[i] for i in range(par_start, par_end)]
                self._params[i][0] *= y_stretch 
            par_start = par_end
    
    def renorm_fits(self):
        par_start = 0
        for i, (num_par, model_type) in enumerate(zip(self._num_params, self._model_types)):
            par_end = par_start + num_par
            if model_type in ['linear']:
                break
            elif model_type in ['polynomial']:
                break
            elif model_type in ['gaussian']:
                break
            elif model_type in ['lorentzian']:
                self._params[i][0] *= func_lorentzian(self._params[i][1], 1, self._params[i][1], self._params[i][2])
            elif model_type in ['voigt']:
                self._params[i][0] *= func_voigt(self._params[i][1], 1, self._params[i][1], self._params[i][2], self._params[i][3])
            elif model_type in ['baskovian']:
                self._params[i][0] *= func_baskovian(self._params[i][1], 1, self._params[i][1], self._params[i][2])
            elif model_type in ['bwf']:
                self._params[i][0] *= func_bwf(self._params[i][1], 1, self._params[i][1], self._params[i][2], self._params[i][3])
                #break
            par_start = par_end
    
    def check_params(self, model_type, params:dict, **kwargs):
        params_list = []
        if model_type in ['linear']:
            keys_required = [key in params for key in ('k', 'd')]
            if all(keys_required):
                params_list = [params[key] for key in ['k', 'd']]
            else:
                print("Linear model needs initial values k and d.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        elif model_type in ['polynomial']:
            for a_i in range(len(params)):
                if 'a'+str(a_i) in params:
                    params_list.append(params['a'+str(a_i)])
                else:
                    print("Polynomial model only supports a_i with i=0,1,2... as initial values.")
        elif model_type in ['gaussian']:
            keys_required = [key in params for key in ('c', 'x0', 'sigma')]
            if all(keys_required):
                params_list = [params[key] for key in ['c', 'x0', 'sigma']]
            else:
                print("Gaussian model needs initial values c, x0 and sigma.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        elif model_type in ['lorentzian']:
            keys_required = [key in params for key in ('c', 'x0', 'gamma')]
            if all(keys_required):
                params_list = [None]*3
                for i_k,key in enumerate(['c', 'x0', 'gamma']):
                    params_list[i_k] = params[key]
                    if key == 'c':
                        params_list[i_k] *= 1/func_lorentzian(params['x0'], 1, params['x0'], params['gamma'])
            else:
                print("Lorentzian model needs initial values c, x0 and gamma.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        elif model_type in ['voigt']:
            keys_required = [key in params for key in ('c', 'x0', 'sigma', 'gamma')]
            if all(keys_required):
                params_list = [None]*4
                for i_k,key in enumerate(['c', 'x0', 'sigma', 'gamma']):
                    params_list[i_k] = params[key]
                    if key == 'c':
                        params_list[i_k] *= 1/func_voigt(params['x0'], 1, params['x0'], params['sigma'], params['gamma'])
            else:
                print("Voigt model needs initial values c, x0, sigma and gamma.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        elif model_type in ['baskovian']:
            keys_required = [key in params for key in ('c', 'x0', 'gamma')]
            if all(keys_required):
                params_list = [None]*3
                for i_k,key in enumerate(['c', 'x0', 'gamma']):
                    params_list[i_k] = params[key]
                    if key == 'c':
                        params_list[i_k] *= 1/func_baskovian(params['x0'], 1, params['x0'], params['gamma'])
            else:
                print("Baskovian model needs initial values c, x0 and gamma.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        elif model_type in ['bwf']:
            keys_required = [key in params for key in ('c', 'x0', 'gamma', 'q')]
            if all(keys_required):
                params_list = [None]*4
                for i_k,key in enumerate(['c', 'x0', 'gamma', 'q']):
                    params_list[i_k] = params[key]
                    if key == 'c':
                        params_list[i_k] *= 1/func_voigt(params['x0'], 1, params['x0'], params['gamma'], params['q'])
            else:
                print("Breit-Wigner Fanon model needs initial values c, x0, gamma and q.")
            if len(params.keys()) > len(keys_required) or keys_required[0] is False:
                print("Too many or wrong keys for this model type.")
        return params_list
    
    def add_fit_func(self, model_type, params:dict, bounds=None):
        self.funcs_update(model_type)
        params_list = self.check_params(model_type, params)
        self.params_update(params_list, bounds)
        self.param_names_update(model_type, params_list)
    def get_params(self, separate:bool=False):
        if separate is False:
            return np.concatenate(self._params)
        else:
            return self._params
    def get_fit_params(self):
        return self.get_params()
    def get_bounds(self):
        return (np.concatenate(self._bounds_low), np.concatenate(self._bounds_high))
    def get_fit_func(self, separate:bool=False):
        if separate is False:
            def func_total(x, *params):
                res = 0.
                par_start = 0 
                for i,func in enumerate(self._funcs):
                    par_end = par_start+self._num_params[i]
                    res += func(x, *params[par_start:par_end])
                    par_start = par_end
                return res
            return func_total
        else:
            return self._funcs
