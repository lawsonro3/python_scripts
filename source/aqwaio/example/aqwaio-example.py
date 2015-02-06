"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import aqwaio as aio
import matplotlib.pyplot as plt
plt.close('all')

# Load AQWA output data file
aq = aio.AqwaOutput(directory='.', outFile='aqwa-example-data.lis')

# Plot diag components of added mass and damping
componentsToPlot = [[0,0],[1,1],[2,2]]
aq.plotAddedMassAndDamping(componentsToPlot)

# Write hydrodynamic data for WEC-Sim
aq.writeWecSimHydroData()

# Write hydrodynamic data to HDF5 file format
aq.writeHdf5()