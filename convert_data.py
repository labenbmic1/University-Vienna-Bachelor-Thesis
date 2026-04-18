"""
@author: Michael Labenbacher
"""
import numpy as np
import scipy.constants

def nm_to_cm1(x, wavelength=633.0):
    x_nm = 10**(7)*(1/wavelength - 1/x) 
    return x_nm

def nm_to_nm_real(wavelength):
    wavelength_ = int(wavelength)
    match wavelength_:
        case 405:
            return 405.
        case 458:
            return 457.935
        case 488:
            return 487.986 
        case 515:
            return 514.532
        case 532: 
            return 530.865
        case 568:
            return 568.188
        case 633:
            return 633.
        case 647:
            return 647.088
        case _:
            return wavelength            

def max_to_one(y, y_max=None, **kwargs):
    # convert max number to 1
    if y_max is None:
        if 'y_max_region' in kwargs:
            y_max = np.max(y[kwargs['y_max_region'][0]:kwargs['y_max_region'][1]+1])
        else:
            y_max = np.max(y)
    y = np.array(y)
    return y/y_max

def nm_to_eV(wavelength):
    return scipy.constants.h*scipy.constants.c*np.power((wavelength * 10**(-9)), -1) * 1/scipy.constants.eV 

def eV_to_nm(energy):
    return scipy.constants.h*scipy.constants.c*np.power((energy * scipy.constants.eV), -1) * 1 / 10**(-9)
