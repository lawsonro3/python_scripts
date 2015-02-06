# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson
"""
import wamitio as wio
import matplotlib.pyplot as plt
plt.close('all')


wamit = wio.WamitOutput(directory='./',outFile='oswec.out')
wamit.plotAddedMassAndDamping(components=[[0,0],[1,1],[2,2]])
wamit.writeHdf5()

