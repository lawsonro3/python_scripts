# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 17:51:21 2014

@author: mlawson
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import collections

class WamitOutput(object):
    '''
    Class to read and interact with WAMIT simulation data
    
    Inputs:
        direcotry: location of the wamit output data
        simName: base name of the wamit simulation files
    Outputs:
        None
    '''
    def __init__(self,directory,outFile):
        self.density = 1000.0
        self.gravity = 9.81 
        self.dir = directory
        self.outFile = outFile


        self.numBodies = 0
        self.bodyNames = {}
        self.bodyPos = {}
        self.period = []
        self.addedMass = {}
        self.addedMassDiag = {}
        self.addedMassAll = {}
        self.radDamping = {}
        self.radDampingDiag = {}
        self.radDampingAll = {}
        self.addedMassAndDampingRaw = {}
        self.files = {}
        self.files['out'] = self.dir + os.path.sep + self.outFile
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
            if "Input from Geometric Data File:" in line:
                self.numBodies = 1
            if " Body number: N=" in line:
                position = []
                self.numBodies += 1
                position.append(self.wamitOutRaw[i+4].split()[2])
                position.append(self.wamitOutRaw[i+4].split()[5])
                position.append(self.wamitOutRaw[i+4].split()[8])
                self.bodyPos[self.numBodies-1] = position
            if "Wave period = infinite" in line:
                self.addedMassInf  = self.wamitOutRaw[i+7:i+7+(6*self.numBodies)**2+1]
            if "Wave period = zero" in line:
                self.addedMassZero = self.wamitOutRaw[i+7:i+7+(6*self.numBodies)**2+1]
            if "Wave period (sec)" in line:
                self.period.append(self.wamitOutRaw[i].split()[4])
                tempFreq = 2.*np.pi/np.float(self.wamitOutRaw[i].split()[4])
                self.addedMassAndDampingRaw[tempFreq] = self.wamitOutRaw[i+7:i+7+(6*self.numBodies)**2]
        self.period = np.array(self.period).astype(float)   
        self.freq = 2.*np.pi/self.period
        self.numFreqs = np.size(self.freq)
        for i,freq in enumerate(self.freq):
            self.addedMassAll[freq]  =  np.array([str(self.addedMassAndDampingRaw[freq][temp]).split()[2] for temp in xrange((6*self.numBodies)**2)]).astype(float).reshape(6*self.numBodies,6*self.numBodies)*self.density
            self.radDampingAll[freq] =  np.array([str(self.addedMassAndDampingRaw[freq][temp]).split()[3] for temp in xrange((6*self.numBodies)**2)]).astype(float).reshape(6*self.numBodies,6*self.numBodies)*self.density*self.freq[i]
        for nb in xrange(self.numBodies):
            addedMass = {}
            radDamping = {}
            addedMassDiag = {}
            radDampingDiag = {}
            for j,freq in enumerate(self.freq):
                addedMass[freq]  = self.addedMassAll[freq][nb*6:nb*6+6]
                radDamping[freq] = self.radDampingAll[freq][nb*6:nb*6+6]
                addedMassDiag[freq] = np.diag(addedMass[freq][:,nb*6:nb*6+6])
                radDampingDiag[freq] = np.diag(radDamping[freq][:,nb*6:nb*6+6])
            self.addedMass[nb] = addedMass
            self.radDamping[nb] = radDamping
            self.addedMassDiag[nb] = addedMassDiag
            self.radDampingDiag[nb] = radDampingDiag
        
    def plotAddedMassAndDamping(self,body=0):
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
        
        am = []
        rad = []
        for i,freq in enumerate(self.freq):
            am.append(self.addedMassDiag[body][freq])
            rad.append(self.radDampingDiag[body][freq])
        am = np.array(am)
        rad = np.array(rad)
        
        for i in xrange(3):
            ax[0].set_title('Diagional Compinent of Added Mass Matrix for Body ' + str(body))
            ax[0].plot(self.freq,am[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
            ax[0].set_ylabel('Added Mass (kg)')
            ax[1].plot(self.freq,rad[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
            ax[1].set_title('Diagional Compinent of Radiation Damping Matrix for Body ' + str(body))
            ax[1].set_xlabel('Wave Frequency (rad/s)')
            ax[1].set_ylabel('Radiation Damping (N-s/m')
            ax[1].legend(loc=0)
            
        plt.show()
#            
class WamitInput(object):
    pass
