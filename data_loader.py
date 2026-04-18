"""
@author: Michael Labenbacher
"""
import numpy as np
path = 'dummy'

def load_data(file_name_add: str="", path:path=path, **kwargs):
    file_name = file_name_add
    with open(path+'/'+file_name, 'r') as file:
        x, y = np.loadtxt(file, delimiter='\t', usecols=(0,1), unpack=True, dtype=np.float64)
    if 'x_invert' in kwargs and kwargs['x_invert'] is True:
        x = np.flip(x)
        y = np.flip(y)
    return np.column_stack((x, y))

def write_data(data, file_name: str="dummy", path:path=path):
    np.savetxt(path+file_name+'.txt', data, delimiter='\t')
        
