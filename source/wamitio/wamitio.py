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
        self.files['out'] = os.path.join(self.dir,self.outFile)
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
                self.bodyNames[self.numBodies-1] = 'body'
            if "Input from Geometric Data Files:" in line:
                for j in xrange(20): # look for bodies within the next 20 lines
                    if "N=" in self.wamitOutRaw[i+j]:
                        self.numBodies += 1
                        self.bodyNames[self.numBodies-1] = self.wamitOutRaw[i+j].split()[-1]
            count = 0 # Counter for bodies
            if " Body number: N=" in line:
                count += 1
                for j in xrange(20): # look for position within the next 20 lines - will only work for wamit files of about 5 bodies
                    if 'XBODY =' in self.wamitOutRaw[i+j]:
                        position = []
                        position.append(self.wamitOutRaw[i+j].split()[2])
                        position.append(self.wamitOutRaw[i+j].split()[5])
                        position.append(self.wamitOutRaw[i+j].split()[8])
                        self.bodyPos[count] = position
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
        
    def plotAddedMassAndDamping(self,bodyToPlot=None):
        '''
        Function to plot the diagional component of added mass and raditation
        dampint
        Inputs:
            None
        Outputs:
            None
        '''
        if bodyToPlot is None:
            bodyToPlot = self.numBodies
        else:
            bodyToPlot = bodyToPlot+1
            
        for body in xrange(bodyToPlot):
            f, ax = plt.subplots(4, sharex=True)
            ax[0].plot()
            ax[1].plot()
            ax[2].plot()
            ax[3].plot()
            
            am = []
            rad = []
            for i,freq in enumerate(self.freq):
                am.append(self.addedMassDiag[body][freq])
                rad.append(self.radDampingDiag[body][freq])
            am = np.array(am)
            rad = np.array(rad)
            
            for i in xrange(3):
                ax[0].set_title('Added Mass for Body ' + str(body))
                ax[0].plot(self.freq,am[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
                ax[0].set_ylabel('Added Mass (kg)')
                ax[1].plot(self.freq,rad[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
                ax[1].set_title('Radiation Damping for Body ' + str(body) + ':' + self.bodyNames[body])
    #            ax[1].set_xlabel('Wave Frequency (rad/s)')
                ax[1].set_ylabel('Radiation Damping (N-s/m)')
                ax[1].legend(loc=0)
            
            for i in xrange(3):
                ax[2].set_title('Added Mass for Body ' + str(body))
                ax[2].plot(self.freq,am[:,i+3],'x-',label='Component (' + str(i+4) + ', ' + str(i+4) + ')')
                ax[2].set_ylabel('Added Mass (kg-m^2)')
                ax[3].plot(self.freq,rad[:,i+3],'x-',label='Component (' + str(i+4) + ':' + str(i+4) + ')')
                ax[3].set_title('Radiation Damping for Body ' + str(body) + ':' + self.bodyNames[body])
                ax[3].set_xlabel('Wave Frequency (rad/s)')
                ax[3].set_ylabel('Radiation Damping (N-m-s/rad)')
                ax[3].legend(loc=0)
                
            plt.show()
#            
class WamitInput(object):
    pass
