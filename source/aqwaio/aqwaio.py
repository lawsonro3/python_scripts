# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import os
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import scipy.io as sio
import re

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
        self.cg = {}
        self.cob = {}
        self.volDisp = {}
        self.periodAllBodies = {}
        self.freqAllBodies = {}
        self.freq = None
        self.minFreq ={}
        self.maxFreq = {}
        self.numFreqs = {}
        self.addedMass = {}
        self.addedMassInfFreq = {}
        self.addedMassZeroFreq = {}
        self.addedMassDiag = {}
        self.addedMassAll = {}
        self.radDamping = {}
        self.radDampingDiag = {}
        self.radDampingAll = {}
        self.addedMassAndDampingRaw = {}
        self.files = {}
        self.files['out'] = self.dir + os.path.sep + self.outFile
        self.files['wecSimHydroData'] = self.dir + os.path.sep + self.outFile[0:-4] + '-wecSimHydroData.mat'
        self.kHeave = {}
        self.kRoll = {}
        self.kPitch = {}
        self.waterPlaneArea = {}
        self.buoyancyForce = {}
        self.kMatrix = {}
        self.exAll = {}
        self.ex = {}
        self.exRe = {}
        self.exIm = {}
        self.exPhase = {}
        self.waterDepth = None
        self.waveDir = {}

        self.readOutFile()

    def readOutFile(self):
        
        with open(self.files['out'],'r') as fid:
            self.outRaw = fid.readlines()
        bodNum = 0
        bodNum2 = 0
        for i, line in enumerate(self.outRaw):
            if 'WATER  DEPTH  . . . . . . . . . . . . . . . . =' in line:
                self.waterDepth = np.array(self.outRaw[i].split())[-1].astype(float)
                self.density = np.array(self.outRaw[i+2].split())[-1].astype(float)
                self.gravity = np.array(self.outRaw[i+4].split())[-1].astype(float)
            if '1. STIFFNESS MATRIX AT THE CENTRE OF GRAVITY' in line:
                cob = []
                self.numBodies += 1
                self.cg[self.numBodies-1] = np.array(self.outRaw[i+2].split())[np.ix_([2,4,6])].astype(float)
                self.volDisp[self.numBodies-1] = np.array(self.outRaw[i+11].split())[-1].astype(float)
                cob.append(np.array(self.outRaw[i+14].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+15].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+16].split())[-1].astype(float))
                self.cob[self.numBodies-1] = np.array(cob)
                self.kHeave[self.numBodies-1]  = np.array(self.outRaw[i+4].split())[-3:].astype(float)
                self.kRoll[self.numBodies-1]   = np.array(self.outRaw[i+5].split())[-3:].astype(float)
                self.kPitch[self.numBodies-1]  = np.array(self.outRaw[i+6].split())[-3:].astype(float) # Fix this!!!
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
                self.periodAllBodies[bodNum] = np.array(period).astype(float)
                self.freqAllBodies[bodNum] = np.array(freq).astype(float)
                self.addedMass[bodNum] = amAll
                self.radDamping[bodNum] = radAll
                self.addedMassDiag[bodNum] = amDiagAll
                self.radDampingDiag[bodNum] = radDiagAll
                self.freq = self.freqAllBodies[0]
                self.maxFreq[bodNum] = np.max(self.freqAllBodies[bodNum])
                self.minFreq[bodNum] = np.max(self.freqAllBodies[bodNum])
                self.numFreqs[bodNum] = np.size(self.freqAllBodies[bodNum])
                self.addedMassInfFreq[bodNum] = self.addedMass[bodNum][self.maxFreq[bodNum]]
                self.addedMassZeroFreq[bodNum] = self.addedMass[bodNum][self.minFreq[bodNum]]
      
                bodNum += 1

            # needs to be fixed to address indexing problems with excitation forces

            if '* * * * H Y D R O D Y N A M I C   P A R A M E T E R S   F O R   S T R U C T U R E'  in line:
                if 'FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY' in self.outRaw[i+4]:
                    self.exAll[bodNum2] = ascii.read(self.outRaw[i+12:i+11+self.numFreqs[0]]) # Change this index from 0 to the correct index
                    temp = self.outRaw[i+11].split()
                    temp2 = float(temp.pop(2))
                    print temp2
                    if temp2 == 0.00:
                        print 'hi there mike'
                        self.waveDir[bodNum2] = temp2
                        self.exAll[bodNum2].add_row(temp)
                        temp = self.exAll[bodNum2].copy()
                        self.exAll[bodNum2][0] = temp[-1]
                        for k,line in enumerate(self.exAll[bodNum2]):
                            if k > 0:
                                self.exAll[bodNum2][k] = temp[k-1]                            
                        self.ex[bodNum2] = np.array([self.exAll[bodNum2].field(2),
                                                   self.exAll[bodNum2].field(4),
                                                   self.exAll[bodNum2].field(6),
                                                   self.exAll[bodNum2].field(8),
                                                   self.exAll[bodNum2].field(10),
                                                   self.exAll[bodNum2].field(12)]).transpose()
                        self.exPhase[bodNum2] = np.array([self.exAll[bodNum2].field(3),
                                                   self.exAll[bodNum2].field(5),
                                                   self.exAll[bodNum2].field(7),
                                                   self.exAll[bodNum2].field(9),
                                                   self.exAll[bodNum2].field(11),
                                                   self.exAll[bodNum2].field(13)]).transpose()
                        self.exRe[bodNum2] = self.ex[bodNum2]*np.cos(np.deg2rad(self.exPhase[bodNum2]))
                        self.exIm[bodNum2]  = self.ex[bodNum2]*np.sin(np.deg2rad(self.exPhase[bodNum2]))
                        bodNum2 += 1

    
    def cutFile1(self):
        self.hydroParmInd = []
        self.addedMassInd = []
        for i, line in enumerate(self.outRaw):
            if '* * * * H Y D R O D Y N A M I C   P A R A M E T E R S   F O R   S T R U C T U R E   1 * * * *' in line:
                self.hydroParmInd.append(i) 
            if 'FROUDE KRYLOV + DIFFRACTION FORCES - VARIATION WITH WAVE DIRECTION' in line:
                self.addedMassInd.append(i)
                
            
                
                
    
    def writeWecSimHydroData(self,bodyNumber=0):
        data = {}
        data['waterDepth'] = self.waterDepth.copy()
        data['waveHeading'] = self.waveDir[bodyNumber]
        data['vol'] = self.volDisp[bodyNumber]
        data['cg'] = self.cg[bodyNumber].copy()
        self.periodAllBodies[bodyNumber].sort()
        data['period'] = self.periodAllBodies[bodyNumber]
        data['linearHyroRestCoef'] = self.kMatrix[bodyNumber]
        data['fAddedMassZero'] = self.addedMassInfFreq[bodyNumber]
        data['fAddedMassInf'] = self.addedMassZeroFreq[bodyNumber]
        amOut = np.zeros((6,6,int(self.numFreqs[bodyNumber])))  
        radOut = np.zeros((6,6,int(self.numFreqs[bodyNumber]))) 
        exReOut = np.zeros((6,int(self.numFreqs[bodyNumber]))) 
        exImOut = np.zeros((6,int(self.numFreqs[bodyNumber])))
        exOut = np.zeros((6,int(self.numFreqs[bodyNumber])))
        exPhaseOut = np.zeros((6,int(self.numFreqs[bodyNumber])))
        keys = self.addedMass[bodyNumber].keys()
        keys.sort()
        for j, freq in enumerate(keys):
            amOut[:,:,j] = self.addedMass[bodyNumber][freq]
            radOut[:,:,j] = self.radDamping[bodyNumber][freq]
        data['fAddedMass'] = amOut
        data['fDamping'] = radOut
        data['fExtRe'] = np.fliplr(self.exRe[bodyNumber].transpose())
        data['fExtIm'] = np.fliplr(self.exIm[bodyNumber].transpose())
        data['fExtMag'] = np.fliplr(self.ex[bodyNumber].transpose())
        data['fExtPhase'] = np.fliplr(self.exPhase[bodyNumber].transpose())
        sio.savemat(self.files['wecSimHydroData'],data)
        
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
            ax[1].set_ylabel('Radiation Damping (N-s/m)')
            ax[1].legend(loc=0)
            
        plt.show()
            
class WamitInput(object):
    pass
