"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson

This is an example of how to use the wamitio module
"""

import wamitio as wio
import matplotlib.pyplot as plt
plt.close('all')
plt.interactive(True)

# Load the data
w = wio.WamitOutput(directory='./',outFile='oswec.out')

# Plot selected components of hydrodynamic coefficinets and excitation force
# Note that python uses zero indexing
comps = [[0,0],[1,1],[2,2]]
w.data[0].plotAddedMassAndDamping(comps)
w.data[1].plotAddedMassAndDamping(comps)
w.data[2].plotAddedMassAndDamping(comps)
w.data[0].plotExcitation([0])

# Save the data in HDF5 format
w.writeHdf5()
