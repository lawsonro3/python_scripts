# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 17:51:21 2014

@author: mlawson
"""
import os
import numpy as np
import hydroData as hd

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

        self.dir = directory
        
        self.density = 1000.

        self.files = {}
        self.files['out'] = os.path.join(self.dir,outFile)
        self.files['hdf5'] = self.files['out'][:-4] + '.h5'
        self.files['wecSim'] = self.files['out'][:-4]
        self.data = {}
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

            wamitOut = fid.readlines()
   
        nBodies = 0 # Total number of bodies
        bodCount = 0 # Counter for bodies
        freqCount = 0
        T = []
        cg = {}
        cb = {}
        name = {}    
        volDisp = {}
        k = {}
        pos = {}
        
        for i, line in enumerate(wamitOut):

            # Read gravity and density
            if 'Gravity:' in line:
                gravity = wamitOut[i].split()[1]
                gravity = np.float(gravity)
            
            if 'Water depth:' in line:
                waterDepth = wamitOut[i].split()[2]
                waterDepth = np.float(waterDepth)

            # If there is one body in the WAMIT run
            if "Input from Geometric Data File:" in line:

                nBodies = 1
                name[0] = 'body'


            
            # If there are two bodies in the WAMIT run
            if "Input from Geometric Data Files:" in line:

                for j in xrange(20): # look for bodies within the next 20 lines

                    if "N=" in wamitOut[i+j]:

                        nBodies += 1
                        name[nBodies-1] = wamitOut[i+j].split()[-1]



            # Read the body positions
            if " Body number: N=" in line:

                for j in xrange(20): # look for position within the next 20 lines - will only work for wamit files of about 5 bodies

                    if 'XBODY =' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        pos[bodCount] = np.array([temp[2],temp[5],temp[8]]).astype(float)
                        
                    if 'Volumes (VOLX,VOLY,VOLZ):' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        volDisp[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)
                        
                    if 'Center of Buoyancy (Xb,Yb,Zb):' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        cb[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)
                        
                    if 'C(3,3),C(3,4),C(3,5):' in wamitOut[i+j]:

                        temp = np.zeros([6,6])
                        temp2 = wamitOut[i+j].split()
                        temp[2,2] = np.float(temp2[1])
                        temp[2,3] = np.float(temp2[2])
                        temp[2,4] = np.float(temp2[3])
                        
                        temp2 = wamitOut[i+j+1].split()
                        temp[3,3] = np.float(temp2[1])
                        temp[3,4] = np.float(temp2[2])
                        temp[3,5] = np.float(temp2[3])
                        
                        temp2 = wamitOut[i+j+2].split()
                        temp[4,4] = np.float(temp2[1])
                        temp[4,5] = np.float(temp2[2])
                        
                        k[bodCount] = temp
                        
                    if 'Center of Gravity  (Xg,Yg,Zg):' in wamitOut[i+j]:
                            
                        temp = wamitOut[i+j].split()
                        cg[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)                                                
                        
                        
                bodCount += 1      


                
            # Inf freq added mass
            if "Wave period = zero" in line:
                
                amInf  = wamitOut[i+7:i+7+(6*nBodies)**2]
                amInf = np.array([amInf[temp].split()[2] for temp in xrange(np.size(amInf))]).astype(float)
                amInf = amInf.reshape(6*nBodies,6*nBodies)


                
            # Zero freq added mass
            if "Wave period = infinite" in line:
                
                amZero = wamitOut[i+7:i+7+(6*nBodies)**2]
                amZero = np.array([amZero[temp].split()[2] for temp in xrange(np.size(amZero))]).astype(float)
                amZero = amZero.reshape(6*nBodies,6*nBodies)

            # Added mass and damping
            if "Wave period (sec)" in line:
                
                T.append(wamitOut[i].split()[4])
                
                am = wamitOut[i+7:i+7+(6*nBodies)**2]
                am = np.array([am[temp].split()[2] for temp in xrange(np.size(am))]).astype(float)
                am = am.reshape(6*nBodies,6*nBodies,1)
                
                rad = wamitOut[i+7:i+7+(6*nBodies)**2]
                rad = np.array([rad[temp].split()[3] for temp in xrange(np.size(rad))]).astype(float)
                rad = rad.reshape(6*nBodies,6*nBodies,1)
                
                ex = wamitOut[i+17+(6*nBodies)**2:i+17+(6*nBodies)**2+6*nBodies]
                ex = np.array([ex[temp].split()[1] for temp in xrange(np.size(ex))]).astype(float)
                ex = ex.reshape(1,6*nBodies)  
                
                phase = wamitOut[i+17+(6*nBodies)**2:i+17+(6*nBodies)**2+6*nBodies]
                phase = np.array([phase[temp].split()[2] for temp in xrange(np.size(phase))]).astype(float)
                phase = phase.reshape(1,6*nBodies)  

                if freqCount is 0:

                    amAll = am
                    radAll = rad
                    exAll = ex
                    phaseAll = phase
                    
                    freqCount = 1

                else:

                    amAll = np.append(amAll,am,axis=2)
                    radAll = np.append(radAll,rad,axis=2)
                    exAll = np.append(exAll,ex,axis=0)
                    phaseAll = np.append(phaseAll,ex,axis=0)
                    
                    
        T = np.array(T).astype(float)
        phaseAll = np.deg2rad(phaseAll)
        exReAll = np.cos(phaseAll)*exAll
        exImAll = np.sin(phaseAll)*exAll
                    
                

        for i in xrange(nBodies):       
            self.data[i] = hd.HydrodynamicData() 
            self.data[i].name = name[i]
            self.data[i].g = gravity
            self.data[i].waterDepth = waterDepth
            self.data[i].rho = self.density            
            self.data[i].nBodies = nBodies
            self.data[i].cg = cg[i]
            self.data[i].cb = cb[i]
            self.data[i].k = k[i]
            self.data[i].pos = pos[i]
            self.data[i].volDisp = volDisp[i]
            
            self.data[i].am.infFreq = amInf[6*i:6+6*i,:]
            self.data[i].am.infFreq = self.data[i].am.infFreq*self.density

            self.data[i].am.zeroFreq = amZero[6*i:6+6*i,:]
            self.data[i].am.zeroFreq = self.data[i].am.zeroFreq*self.density
            
            self.data[i].T = T
            self.data[i].w = 2.0*np.pi/self.data[i].T

            self.data[i].am.all = amAll[6*i:6+6*i,:,:]
            self.data[i].am.all = self.data[i].am.all*self.density
            
            self.data[i].rd.all = radAll[6*i:6+6*i,:,:]
            for j in xrange(np.shape(self.data[i].rd.all)[2]):
                self.data[i].rd.all[:,:,j] = self.data[i].rd.all[:,:,j]*self.density*self.data[i].w[j]
                
            self.data[i].ex.mag = exAll[:,6*i:6+6*i]
            self.data[i].ex.phase = phaseAll[:,6*i:6+6*i]
            self.data[i].ex.re = exReAll[:,6*i:6+6*i]
            self.data[i].ex.im = exImAll[:,6*i:6+6*i]



    def writeWecSimHydroData(self):
        hd.writeWecSimHydroData(self.data,self.files['wecSim'])    
        
        
        
    def writeHdf5(self):
        hd.writeHdf5(self.data,self.files['hdf5'])
       
       
       
    def plotAddedMassAndDamping(self,components):
        hd.plotAddedMassAndDamping(self.data,components)

