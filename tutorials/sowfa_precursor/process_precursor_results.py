#! /Users/mlawson/env/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import imp
import os

# import precursor functions
import python_scripts.openfoam.sowfa_precursor as sowfa
from python_scripts.openfoam.read import read_input as read_set_up
from python_scripts.materials_properties.physical_constants import Constants as const

plt.close('all')

setUp = read_set_up('/Users/mlawson/scratch/setUp')# Read data from SOWFA 'setUp' file
inputData = setUp.data
inputData['avg_time'] = 18000                      # list the times in the averaging directory
inputData['avg_width'] = 2000                      # time to average about (s)
inputData['mov_avg_width'] = 2000                  # averaging width (s)
inputData['times'] = '0'
imputData = {}
inputData['base_dir'] = '/Users/mlawson/scratch'
inputData['avg_dir'] = 'postProcessing/averaging'
inputData['time_dir'] = os.path.join(inputData['base_dir'], inputData['avg_dir'], inputData['times'])
inputData['zLevel'] = 90.0                         # height to plot average velocity vs. time (m)
inputData['alpha'] = 0.17                          # I think this can be deleted
inputData['heights'] = [0.0, inputData['windHeight']/2, inputData['windHeight'], 3*inputData['windHeight']/2, 2*inputData['windHeight'], 3*inputData['windHeight'], 4*inputData['windHeight'], 5*inputData['windHeight']]


# ?? Whart is this ??
#uStarAvg = sowfa.uStar_avg_cell(baseDir,avg_dir,time,avg_time,avg_width)
inputData['uStarAvg'] = 0.72

inputData['zCell'] = sowfa.read_h_levels_cell(inputData)
inputData['ziAvg'] = sowfa.theta_w_avg_cell(inputData)

inputData = sowfa.Umean_avg_nonnormalized(inputData)
inputData['Tmag'] = sowfa.Tmean_avg_nonnormalized(inputData)
inputData = sowfa.variances_avg_cell(inputData)
sowfa.Umean_h(inputData)


# plot other stuff
z = np.linspace(inputData['zMin'],inputData['zMax'],1000)
U = inputData['U0Mag']*(z/inputData['windHeight'])**inputData['alpha']

plt.figure(num='U vs. z')
plt.plot(U,z,'r-',label='U')
plt.plot([2,12],[inputData['windHeight'],inputData['windHeight']],'k--')


#wStarAvg = ((g/inputData['TRef'])*Q_s*inputData['ziAvg'])**(1./3.) # Jen's origional with Qs=0
wStarAvg = ((const().g/inputData['TRef'])*inputData['heatingRate']*inputData['ziAvg'])**(1./3.)
#L = -(inputData['TRef']*inputData['uStarAvg']**3)/(g*Q_s*0.4)
#neg_zi_L = -inputData['ziAvg']/L
tau_uStar = inputData['ziAvg']/inputData['uStarAvg']
tau_wStar = inputData['ziAvg']/wStarAvg

U_top = np.sqrt( inputData['Uvec'][3,0]**2 + inputData['Uvec'][3,1]**2 )
U_hub = np.sqrt( inputData['Uvec'][2,0]**2 + inputData['Uvec'][2,1]**2 )
U_bot = np.sqrt( inputData['Uvec'][1,0]**2 + inputData['Uvec'][1,1]**2 )

U_shear = (U_top - U_bot)/(inputData['heights'][2] - inputData['heights'][0])
dir_shear = (inputData['dir'][2] - inputData['dir'][0])/(inputData['heights'][2]-inputData['heights'][0])

print('Top = ', U_top)
print('Hub = ', U_hub)
print('Bot = ', U_bot)

print(inputData)
