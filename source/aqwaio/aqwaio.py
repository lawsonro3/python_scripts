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
    def __init__(self,directory,outFile):
        self.density = 1000.0
        self.gravity = 9.81 
        self.dir = directory
        self.outFile = outFile


        self.numBodies = 0
        self.bodyNames = {}
        self.cog = {}
        self.cob = {}
        self.volDisp = {}
        self.period = {}
        self.freq = {}
        self.addedMass = {}
        self.addedMassDiag = {}
        self.addedMassAll = {}
        self.radDamping = {}
        self.radDampingDiag = {}
        self.radDampingAll = {}
        self.addedMassAndDampingRaw = {}
        self.files = {}
        self.files['out'] = self.dir + os.path.sep + self.outFile
        self.kHeave = {}
        self.kRoll = {}
        self.kPitch = {}
        self.waterPlaneArea = {}
        self.buoyancyForce = {}
        self.kMatrix = {}

        self.readOutFile()

        
        
    def readOutFile(self):
        
        with open(self.files['out'],'r') as fid:
            self.outRaw = fid.readlines()
        cob = []
        bodNum = 0
        for i, line in enumerate(self.outRaw):
            if '1. STIFFNESS MATRIX AT THE CENTRE OF GRAVITY' in line:
                self.numBodies += 1
                self.cog[self.numBodies-1] = np.array(self.outRaw[i+2].split())[np.ix_([2,4,6])].astype(float)
                self.volDisp[self.numBodies-1] = np.array(self.outRaw[i+11].split())[-1].astype(float)
                cob.append(np.array(self.outRaw[i+14].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+15].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+16].split())[-1].astype(float))
                self.cob[self.numBodies-1] = np.array(cob)
                self.kHeave[self.numBodies-1] = np.array(self.outRaw[i+6].split())[-3:].astype(float)
                self.kRoll[self.numBodies-1] = np.array(self.outRaw[i+7].split())[-3:].astype(float)
#                self.kPitch[self.numBodies-1] = np.array(self.outRaw[i+8].split())[-3:].astype(float) # Fix this!!!
                self.waterPlaneArea = np.array(self.outRaw[i+19].split())[-1].astype(float)
            
            if 'AT THE FREE-FLOATING EQUILIBRIUM POSITION' in line:
                runLoop = True
                period = []
                freq = []
                amAll = {}
                radAll = {}
                kMatrix = []
                amDiagAll = {}
                radDiagAll = {}
                
            
                ind = 31
                
                self.buoyancyForce[bodNum] = float(self.outRaw[i+3].split()[-1])
                
                kMatrix.append(np.array(self.outRaw[i+16].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+18].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+20].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+22].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+24].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+26].split()[1:]).astype(float))
                
                self.kMatrix[bodNum] = np.array(kMatrix)
                
                while runLoop is True:  
                    am = []
                    rad = []
                    if self.outRaw[i+ind+45].split()[0] != 'WAVE':
                        runLoop = False
                        
                    period.append(self.outRaw[i+ind].split()[3])
                    freq.append(self.outRaw[i+ind].split()[-1])
                    
                    am.append(np.array(self.outRaw[i+ind+10].split()[1:]).astype(float))
                    am.append(np.array(self.outRaw[i+ind+12].split()[1:]).astype(float))
                    am.append(np.array(self.outRaw[i+ind+14].split()[1:]).astype(float))
                    am.append(np.array(self.outRaw[i+ind+16].split()[1:]).astype(float))
                    am.append(np.array(self.outRaw[i+ind+18].split()[1:]).astype(float))
                    am.append(np.array(self.outRaw[i+ind+20].split()[1:]).astype(float))
                    
                    
                    rad.append(np.array(self.outRaw[i+ind+30].split()[1:]).astype(float))
                    rad.append(np.array(self.outRaw[i+ind+32].split()[1:]).astype(float))
                    rad.append(np.array(self.outRaw[i+ind+34].split()[1:]).astype(float))
                    rad.append(np.array(self.outRaw[i+ind+36].split()[1:]).astype(float))
                    rad.append(np.array(self.outRaw[i+ind+38].split()[1:]).astype(float))
                    rad.append(np.array(self.outRaw[i+ind+40].split()[1:]).astype(float))

                    amAll[float(freq[-1])] = np.array(am)
                    radAll[float(freq[-1])] = np.array(rad)
                    amDiagAll[float(freq[-1])] = np.diag(np.array(am))
                    radDiagAll[float(freq[-1])] = np.diag(np.array(rad))
                    
                    ind += 45
                amAll[freq[-1]] = am
                self.period[bodNum] = np.array(period).astype(float)
                self.freq[bodNum] = np.array(freq).astype(float)
                self.addedMass[bodNum] = amAll
                self.radDamping[bodNum] = radAll
                self.addedMassDiag[bodNum] = amDiagAll
                self.radDampingDiag[bodNum] = radDiagAll
                
                bodNum += 1
                    
                    
                     
        
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
        for i,freq in enumerate(self.freq[body]):
            am.append(self.addedMassDiag[body][freq])
            rad.append(self.radDampingDiag[body][freq])
        am = np.array(am)
        rad = np.array(rad)
        
        for i in xrange(3):
            ax[0].set_title('Diagional Compinent of Added Mass Matrix for Body ' + str(body))
            ax[0].plot(self.freq[body],am[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
            ax[0].set_ylabel('Added Mass (kg)')
            ax[1].plot(self.freq[body],rad[:,i],'x-',label='Component (' + str(i+1) + ', ' + str(i+1) + ')')
            ax[1].set_title('Diagional Compinent of Radiation Damping Matrix for Body ' + str(body))
            ax[1].set_xlabel('Wave Frequency (rad/s)')
            ax[1].set_ylabel('Radiation Damping (N-s/m)')
            ax[1].legend(loc=0)
            
        plt.show()
            
class WamitInput(object):
    pass
