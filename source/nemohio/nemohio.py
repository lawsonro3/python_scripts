# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 15:28:52 2014

@author: mlawson

This module reads data from Nemoh simulations
"""
import numpy as np
import matplotlib.pyplot as plt

class NemohOutput(object):
    def __init__(self,outputDir,plotData=False):
        self.w = []
        self.addedMass = {}
        self.radiationDamping = {}
        self.outputDir = outputDir
        
        self.readCMCA()
        
        if plotData is True:
            self.plotAddedMassAndDamping()
        
    def readCMCA(self):
        '''
        Function to read CM.dat created by Nemoh
        '''        
        # load  data files
        with open(self.outputDir + '/results/CM.dat') as fid :
            linesCM = fid.readlines()        
        with open(self.outputDir + '/results/CA.dat') as fid :
            linesCA = fid.readlines()

        # Read the number of frequencies
        self.numFreqs = int(linesCM[0].split()[-1])

        # Read the Frequencies, the added mass matrix, and the radiation damping matrix at each frequency
        for i in xrange(self.numFreqs):
            self.w.append(float(linesCM[1+i*7].replace('\n','')))
            self.addedMass[self.w[i]] = [temp.replace('\n','') for temp in linesCM[2+i*7:8+i*7]]
            self.addedMass[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.addedMass[self.w[i]]])        
            self.radiationDamping[self.w[i]] = [temp.replace('\n','') for temp in linesCA[2+i*7:8+i*7]]
            self.radiationDamping[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.radiationDamping[self.w[i]]])
            
        self.addedMassDiag = np.array([np.diag(temp) for temp in self.addedMass.values()])
        self.radiationDampingDiag = np.array([np.diag(temp) for temp in self.radiationDamping.values()])
            
    def plotAddedMassAndDamping(self):
        plt.figure('AddedMass')
        plt.figure('RadiationDamping')        
        
        for i in xrange(6):
            plt.figure('AddedMass')
            plt.plot(self.w,self.addedMassDiag[:,i],label='Component (' + str(i) + ', ' + str(i) + ')')
            plt.title('Diagional Compinent of Added Mass Matrix')
            plt.xlabel('Wave Frequency (rad)')
            plt.ylabel('Added Mass (kg)')
            plt.legend()
            
            plt.figure('RadiationDamping') 
            plt.plot(self.w,self.radiationDampingDiag[:,i],label='Component (' + str(i) + ', ' + str(i) + ')')
            plt.title('Diagional Compinent of Radiation Damping Matrix')
            plt.xlabel('Wave Frequency (rad)')
            plt.ylabel('Radiation Damping')
            plt.legend()
            

        
            
test = NemohOutput('/Users/mlawson/Applications/nemoh/matlabRoutines/simulationTest',plotData=True)


        
