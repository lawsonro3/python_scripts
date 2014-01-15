# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:33:01 2012

@author: mlawson
"""
from pylab import *

def BackwardsEuler(data):
    '''
    Compute the first order backwards derivative of a dataset
    '''   
    d_dx = zeros(len(data))
    for i in range(len(data)):
        if i == 0 or i == len(data)-1:
            d_dx[i] = 0 # assume the slope is 0 at the boundary
        else:
            d_dx[i] = data[i]-data[i-1]
    return d_dx