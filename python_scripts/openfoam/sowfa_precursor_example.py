import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# import imp
import os
from glob import glob

# User functions
import read
import sowfa_precursor

inputData =  sowfa_precursor.create(dir='/Users/mlawson/Peregrine/windsim/wake_steering/stableABLRuns/infPer_0.001m_5m',time_dir=25000)
inputData['avg_time'] = 29000                      # time to average about (s)
inputData['avg_width'] = 2000                      # average width (s)
inputData['zLevel'] = 90.0                         # height to plot average velocity vs. time (m)
inputData['uStarAvg'] = 0.72

setUp = read.read_input(inputData['setUpFile'])
inputData.update(setUp.data)

inputData['heights'] = np.array([0,1/2,1,3/2,2,3,4,5])*inputData['windHeight']


inputData['zCell'] = np.array(open(inputData['hLevelCellFile'],'r').read().split()).astype(np.float)

# print(inputData)


inputData = sowfa_precursor.theta_w_avg_cell(inputData)

inputData = sowfa_precursor.Umean_avg_nonnormalized(inputData)

inputData = sowfa_precursor.Tmean_avg_nonnormalized(inputData)

inputData = sowfa_precursor.variances_avg_cell(inputData)





