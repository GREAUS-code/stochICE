# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:57:02 2024

@author: dugj2403
"""

import os

path=os.getcwd()


files=['CD1TEST.txt',
       'CD1TEST_intp.txt',
       'CD1TEST_no_intp.txt',
       'DOUT1.TXT',
       'DOUT2.TXT',
       'DOUT7.TXT',
       'INPT10.TXT',
       'INPT11.TXT',
       'INPT101.TXT',
       'INPT111.TXT',
       'psplot.ps',
       'TAPE5.txt', 'TAPE6.txt','TAPE10.txt',
       'TAPE12.txt','TAPE13.txt','TAPE14.txt',
       'TAPE15.txt','TAPE16.txt','TAPE17.txt',
       'TAPE66.txt','TAPE70.txt','TAPE71.txt',
       'TAPE72.txt','TAPE73.txt','TAPE5.txt',
       'intpxs.txt','intpxs1.txt',
       'TESTCD2.txt']

for file in files:

    try:
        os.remove(path +'\\' + file)
    except FileNotFoundError:
        pass
        




