# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 17:51:21 2014

@author: mlawson
"""
import os
import numpy as np
import matplotlib.pyplot as plt

class WamitOutput(object):
    '''
    Class to read and interact with WAMIT simulation data
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
        '''
        Function to read WAMIT output file into the class
        Inputs:
            None
        Outputs:
            None
        '''
        with open(self.files['out'],'r') as fid:
            self.wamitOutRaw = fid.readlines()
        
        for i, line in enumerate(self.wamitOutRaw):
            if "Wave period = infinite" in line:
                self.addedMassInf = self.wamitOutRaw[i+7:i+7+self.numBodies*37]
            if "Wave period = zero" in line:
                self.addedMassZero = self.wamitOutRaw[i+7:i+7+self.numBodies*37]
            if "Wave period (sec)" in line:
                self.period.append(self.wamitOutRaw[i].split()[4])
                tempFreq = 2.*np.pi/np.float(self.wamitOutRaw[i].split()[4])
                self.addedMassAndDampingRaw[tempFreq] = self.wamitOutRaw[i+7:i+7+self.numBodies*36]
        self.period = np.array(self.period).astype(float)   
        self.freq = 2.*np.pi/self.period
        for i,freq in enumerate(self.freq):
            self.addedMass[freq]  =  np.array([str(self.addedMassAndDampingRaw[freq][temp]).split()[2] for temp in xrange(self.numBodies*36)]).astype(float).reshape([6,6])*self.density
            self.radDamping[freq] =  np.array([str(self.addedMassAndDampingRaw[freq][temp]).split()[3] for temp in xrange(self.numBodies*36)]).astype(float).reshape([6,6])*self.density
            self.radDamping[freq] = self.radDamping[freq]*self.freq[i]
            
            
        addedMassDiag = []
        radiationDampingDiag = []
        for i,freq in enumerate(self.freq):
            addedMassDiag.append(np.diag(self.addedMass[freq]))
            radiationDampingDiag.append(np.diag(self.radDamping[freq]))
            
        self.addedMassDiag = np.array(addedMassDiag)
        self.radiationDampingDiag = np.array(radiationDampingDiag)
        
    def plotAddedMassAndDamping(self):
        '''
        Function to plot the diagional component of added mass and raditation
        dampint
        Inputs:
            None
        Outputs:
            None
        '''
        f, ax = plt.subplots(2, sharex=True)
        ax[0].plot()
        ax[1].plot()
        
        for i in xrange(3):
            ax[0].set_title('Diagional Compinent of Added Mass Matrix')
            ax[0].plot(self.freq,self.addedMassDiag[:,i],'x-',label='Component (' + str(i) + ', ' + str(i) + ')')
            ax[0].set_ylabel('Added Mass (kg)')
            ax[1].plot(self.freq,self.radiationDampingDiag[:,i],'x-',label='Component (' + str(i) + ', ' + str(i) + ')')
            ax[1].set_title('Diagional Compinent of Radiation Damping Matrix')
            ax[1].set_xlabel('Wave Frequency (rad/s)')
            ax[1].set_ylabel('Radiation Damping')
            ax[1].legend(loc=0)
            
        plt.show()
            
class WamitInput(object):
    pass
