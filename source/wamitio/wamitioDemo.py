# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson
"""
import wamitio as wio

wamit = wio.WamitOutput(directory='/Users/mlawson/Applications/nemoh/matlabRoutines/nonsymmetrical-wamit',simName='skewed-shpere')
wamit.plotAddedMassAndDamping()