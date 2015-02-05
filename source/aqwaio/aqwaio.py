# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import os
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import hydroData as hd

class AqwaOutput(object):
    '''
    Class to read and interact with AQWA simulation data
    Inputs:
        direcotry: location of the wamit output data
        simName: base name of the wamit simulation files
    Outputs:
        None
    '''
    def __init__(self,directory,outFile ):

        self.dir = directory
        
        self.files = {}
        self.files['out'] = os.path.join(self.dir,outFile)
        self.files['hdf5'] = os.path.join(self.dir,outFile[0:-4] + '.h5')
        self.files['wecSim'] = os.path.join(self.dir,outFile[0:-4])
        
        self.data = {}        
        
        self.nBodies = 0
        self.gravity = None
        self.density = None
        self.waterDepth = None
 
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
                self.data[self.nBodies] = hd.HydrodynamicData()
                cob = []                
                self.data[self.nBodies].cg = np.array(self.outRaw[i+2].split())[np.ix_([2,4,6])].astype(float)
                self.data[self.nBodies].volDisp = np.array(self.outRaw[i+11].split())[-1].astype(float)
                cob.append(np.array(self.outRaw[i+14].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+15].split())[-1].astype(float))
                cob.append(np.array(self.outRaw[i+16].split())[-1].astype(float))
                self.data[self.nBodies].cb = np.array(cob)
                self.data[self.nBodies].wpArea = np.array(self.outRaw[i+19].split())[-1].astype(float)
                self.data[self.nBodies].waterDepth = self.waterDepth
                self.data[self.nBodies].gravity = self.gravity 
                self.data[self.nBodies].density = self.gravity 
                self.nBodies += 1
            
            if 'AT THE FREE-FLOATING EQUILIBRIUM POSITION' in line:
                runLoop = True
                period = []
                freq = []
                kMatrix = []
                
                ind = 31
                
                self.data[bodNum].buoyForce = float(self.outRaw[i+3].split()[-1])
                
                kMatrix.append(np.array(self.outRaw[i+16].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+18].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+20].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+22].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+24].split()[1:]).astype(float))
                kMatrix.append(np.array(self.outRaw[i+26].split()[1:]).astype(float))
                
                self.data[bodNum].k = np.array(kMatrix)                
                count = 0
                
                while runLoop is True:  
                    am = []
                    rd = []
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
                    am = np.array(am).reshape(6,6,1)
                    
                    rd.append(np.array(self.outRaw[i+ind+30].split()[1:]).astype(float))
                    rd.append(np.array(self.outRaw[i+ind+32].split()[1:]).astype(float))
                    rd.append(np.array(self.outRaw[i+ind+34].split()[1:]).astype(float))
                    rd.append(np.array(self.outRaw[i+ind+36].split()[1:]).astype(float))
                    rd.append(np.array(self.outRaw[i+ind+38].split()[1:]).astype(float))
                    rd.append(np.array(self.outRaw[i+ind+40].split()[1:]).astype(float))
                    rd = np.array(rd).reshape(6,6,1)

                    if count is 0:
                        amAll = am
                        rdAll = rd
                    else:
                        amAll = np.append(amAll,am,axis=2)
                        rdAll = np.append(rdAll,rd,axis=2)
                    
                    ind += 45
                    count += 1
                    
                self.data[bodNum].T = np.array(period).astype(float)
                self.data[bodNum].w = np.array(freq).astype(float)
                self.data[bodNum].am.all = amAll
                self.data[bodNum].rd.all = rdAll
                self.data[bodNum].am.infFreq = self.data[bodNum].am.all[:,:,-1]
                self.data[bodNum].am.zeroFreq = self.data[bodNum].am.all[:,:,0]

                self.data[bodNum].nW = np.size(self.data[bodNum].w)

                bodNum += 1

            if '* * * * H Y D R O D Y N A M I C   P A R A M E T E R S   F O R   S T R U C T U R E'  in line:
                if 'FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY' in self.outRaw[i+4]:
                    temp3 = ascii.read(self.outRaw[i+12:i+11+np.size(self.data[0].w)]) # Change this index from 0 to the correct index
                    temp = self.outRaw[i+11].split()
                    temp2 = float(temp.pop(2))
                    exMag = []
                    exPhase = []
                    exAll = {}
                    if temp2 == 180:
                        exAll[bodNum2] = temp3
                        exAll[bodNum2].add_row(temp)
                        temp = exAll[bodNum2].copy()
                        exAll[bodNum2][0] = temp[-1]
                        for k,line in enumerate(exAll[bodNum2]):
                            if k > 0:
                                exAll[bodNum2][k] = temp[k-1]  
                        for m,freq in enumerate(exAll[bodNum2].field(1)):
                            exMag.append(np.array([exAll[bodNum2].field(2)[m],
                                                       exAll[bodNum2].field(4)[m],
                                                       exAll[bodNum2].field(6)[m],
                                                       exAll[bodNum2].field(8)[m],
                                                       exAll[bodNum2].field(10)[m],
                                                       exAll[bodNum2].field(12)[m]]))
                            exPhase.append(np.array([exAll[bodNum2].field(3)[m],
                                                       exAll[bodNum2].field(5)[m],
                                                       exAll[bodNum2].field(7)[m],
                                                       exAll[bodNum2].field(9)[m],
                                                       exAll[bodNum2].field(11)[m],
                                                       exAll[bodNum2].field(13)[m]]))

                        self.data[bodNum2].waveDir = temp2                               
                        self.data[bodNum2].ex.mag = np.array(exMag)
                        self.data[bodNum2].ex.phase = np.array(exPhase)
                        self.data[bodNum2].ex.re = self.data[bodNum2].ex.mag*np.cos(np.deg2rad(self.data[bodNum2].ex.phase))
                        self.data[bodNum2].ex.im  = self.data[bodNum2].ex.mag*np.sin(np.deg2rad(self.data[bodNum2].ex.phase))
                        bodNum2 += 1
        
    
    def writeWecSimHydroData(self):
        hd.writeWecSimHydroData(self.data,self.files['wecSim'])     
        
    def writeHdf5(self):
        hd.writeHdf5(self.data,self.files['hdf5'])
        
    def plotAddedMassAndDamping(self,components):
        hd.plotAddedMassAndDamping(self.data,components)
