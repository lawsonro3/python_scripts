"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson

This is an example of how to use the wamitio module
"""

import wamitio as wio
import matplotlib.pyplot as plt
plt.close('all')

# Load the data
w = wio.WamitOutput(directory='./',outFile='oswec.out')

# Plot the 1,1, 2,2 and 3,3 components of the added mass matrix
# Note that python uses zero indexing
componentsToPlot = [[0,0],[1,1],[2,2]]
w.plotAddedMassAndDamping(components=componentsToPlot)

# Save the data for use in WEC-Sim
w.writeWecSimHydroData()

# Save the data in HDF5 format
w.writeHdf5()
