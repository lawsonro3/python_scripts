# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import os
import numpy as np
import matplotlib.pyplot as plt

class AqwaOutput(object):
    '''
    Class to read and interact with AQWA simulation data
    Inputs:
        direcotry: location of the wamit output data
        simName: base name of the wamit simulation files
    Outputs:
        None
    '''
    def __init__(self,directory,simName):
        self.density = 1000.0
        self.gravity = 9.81 
        self.dir = directory
        self.simName = simName
        self.numBodies = 1
        
        self.period = []
        self.addedMass = {}
        self.radDamping = {}
        self.addedMassAndDampingRaw = {}
        self.files = {}
        self.files['out'] = self.dir + os.path.sep + self.simName + '.out'
        self.readOutFile()
        
    def readOutFile(self):
        pass # need to write
  
    def plotAddedMassAndDamping(self):
        pass # need to write
            
class WamitInput(object):
    pass
