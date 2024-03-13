# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import h5py
import pandas as pd
import numpy as np
filename = r'C:\Users\Jason\Dropbox\PythonTools\stochICE\Simple_RIVICE_debugging\RIVICE_debug.p01.hdf'

# reread = pd.read_hdf(filename) 


with h5py.File(filename, "r+") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    # print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    # a_group_key = list(f.keys())[3]

    # # get the object type for a_group_key: usually group or dataset
    # print(type(f[a_group_key])) 

    # # If a_group_key is a group name, 
    # # this gets the object names in the group and returns as a list
    # data = list(f[a_group_key])

    # # If a_group_key is a dataset name, 
    # # this gets the dataset values and returns as a list
    # data = list(f[a_group_key])
    # # preferred methods to get dataset values:
    # bob=f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()]
    # ds_obj = f[a_group_key]      # returns as a h5py dataset object
    # ds_arr = f[a_group_key][()]  # returns as a numpy array

    newWSE = np.asarray([69.574814, 69.53983, 69.51734, 69.49125, 69.47724, 69.45276, 69.44871, 69.47081, 69.45508, 69.43802, 69.41144, 69.38965, 69.38325, 69.39193], dtype='float32')
    
    bob=f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()]
    f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()]=newWSE
    billy=f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()]
    # print(f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()])
    f.close()


# f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()][0]
newWSE
bob[()]
