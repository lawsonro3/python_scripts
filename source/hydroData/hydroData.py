# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:50:35 2015

@author: mlawson
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
    
class HydrodynamicCoefficients(object):
    def __init__(self):
        self.all        = np.array([])
        self.infFreq    = np.array([])
        self.zeroFreq   = np.array([])

class HydrodynamicExcitation(object):
    def __init__(self):
        self.re         = np.array([])
        self.im         = np.array([])
        self.mag        = np.array([])
        self.phase      = np.array([])
        
class HydrodynamicData(object):
    def __init__(self):
        self.rho = None
        self.g = None      
        self.files = {}
        self.nBodies = 0                    # Number of bodies in simulation
        self.cg = {}                        # Center of gravity
        self.cb = {}                        # Center of buoyancy
        self.volDisp = {}                   # Volume displacement
        self.T = {}                         # Wave period
        self.w = {}                         # Wave freq
        self.am = HydrodynamicCoefficients()# Added mass 
        self.rd = HydrodynamicCoefficients()# Radiation damping
        self.wpArea = {}                    # Water plane area          
        self.buoyForce = {}                 # Buoyanch force at equelibrium
        self.k = {}                         # Hydrostatic stifness matrix
        self.ex = HydrodynamicExcitation()  # Excitation coeffs
        self.waterDepth = None              # Water depth
        self.waveDir = {}                   # Wave direction
        self.name = None                    # Name of the body
        
def writeHdf5(data,outFile):
        try:
            import h5py
        except:
            raise Exception('The h5py module must be installed to used the writeHdf5 functionality.')
            
        with h5py.File(outFile, "w") as f:        
            for i in range(np.size(data.keys())):
                
                per = f.create_dataset('body' + str(i) + '/period',data=data[i].T)
                per.attrs['units'] = 's'                
                
                freq = f.create_dataset('body' + str(i) + '/frequency',data=data[i].w)
                freq.attrs['units'] = 'rad/s'                
                
                kMat = f.create_dataset('body' + str(i) + '/stiffnessMatrix',data=data[i].k)
                kMat.attrs['units'] = ''
                
                cg = f.create_dataset('body' + str(i) + '/centerOfGravity',data=data[i].cg)
                cg.attrs['units'] = 'm'

                cb = f.create_dataset('body' + str(i) + '/centerOfBuoyancy',data=data[i].cb)
                cb.attrs['units'] = 'm'
                
                exMag = f.create_dataset('body' + str(i) + '/excitationForce/magnitude',data=data[i].ex.mag)
                exMag.attrs['units'] = ''
                
                exPhase = f.create_dataset('body' + str(i) + '/excitationForce/phase',data=data[i].ex.phase)
                exPhase.attrs['units'] = 'rad'
                
                exRe = f.create_dataset('body' + str(i) + '/excitationForce/imaginaryComponent',data=data[i].ex.im)
                exRe.attrs['units'] = ''

                exIm = f.create_dataset('body' + str(i) + '/excitationForce/realComponent',data=data[i].ex.re)
                exIm.attrs['units'] = ''

                # Write added mass information                
                amInf = f.create_dataset('body' + str(i) + '/addedMassCoefficients/infFreq',data=data[i].am.infFreq)
                amInf.attrs['units for translational degrees of freedom'] = 'kg'
                
                am = f.create_dataset('body' + str(i) + '/addedMassCoefficients/discreteFeqs',data=data[i].am.all)
                am.attrs['units for translational degrees of freedom'] = 'kg'                
                am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
                
                for m in xrange(np.shape(data[i].am.all)[0]):
                    for n in xrange(np.shape(data[i].am.all)[1]):
                        amComp = f.create_dataset('body' + str(i) + '/addedMassCoefficients/components/' + str(m) + ',' + str(n),data=data[i].am.all[m,n,:])
                        amComp.attrs['units'] = ''

                        radComp = f.create_dataset('body' + str(i) + '/radiationDampingCoefficients/components/' + str(m) + ',' + str(n),data=data[i].rd.all[m,n,:])
                        radComp.attrs['units'] = ''
                
                rad = f.create_dataset('body' + str(i) + '/radiationDampingCoefficients/discreteFeqs',data=data[i].rd.all)
                rad.attrs['units'] = ''
                
                wDepth = f.create_dataset('body' + str(i) + '/waterDepth',data=data[i].waterDepth)
                wDepth.attrs['units'] = 'm'
                
                try:
                    waveHead = f.create_dataset('body' + str(i) + '/waveDirection',data=np.deg2rad(data[i].waveDir))
                    waveHead.attrs['units'] = 'rad'
                except:
                    pass
                
                
                vol = f.create_dataset('body' + str(i) + '/displacedVolume',data=data[i].volDisp)
                vol.attrs['units'] = 'm^3'

                
                g = f.create_dataset('body' + str(i) + '/gravity',data=data[i].g)
                g.attrs['units'] = 'm/s^2'
    
def writeWecSimHydroData(data,outFile):
        for i in range(np.size(data.keys())):
            curData = data[i]
            out = {}
            out['waterDepth'] = curData.waterDepth
            out['waveHeading'] = curData.waveDir
            out['vol'] = curData.volDisp
            out['cg'] = curData.cg
            out['period'] = curData.T[::-1]
            out['linearHyroRestCoef'] = curData.k
            out['fAddedMassZero'] = curData.am.infFreq
            out['fAddedMass'] = curData.am.all[:,:,::-1]
            out['fDamping'] = curData.am.all[:,:,::-1]
            out['fExtRe'] = curData.ex.re[::-1,:].transpose()
            out['fExtIm'] = curData.ex.im[::-1,:].transpose()
            out['fExtMag'] = curData.ex.mag[::-1,:].transpose()
            out['fExtPhase'] = curData.ex.phase[::-1,:].transpose()
            
            outFileName = outFile+'-body' + str(i) +'.mat'
            sio.savemat(outFileName,out)
            
def plotAddedMassAndDamping(data,components):
    '''
    Function to plot the diagional component of added mass and raditation
    damping - this could be singificantly improved
    Inputs:
        A list of components to plot
    Outputs:
        None
    '''
    for body in xrange(np.size(data.keys())):
        
        f, ax = plt.subplots(2, sharex=True)
        ax[0].plot()
        ax[0].set_title('Added mass for body ' + str(body) + ': ' + str(data[body].name))    
        ax[0].set_ylabel('Added Mass')
        
        ax[1].plot()
        ax[1].set_title('Radiation damping body ' + str(body) + ': ' + str(data[body].name))
        ax[1].set_xlabel('Wave Frequency (rad/s)')
        ax[1].set_ylabel('Radiation Damping')
        
        for i,comp in enumerate(components):
            
            x = comp[0]
            y = comp[1]
            w = data[body].w
            rd = data[body].rd.all[x,y,:]
            am = data[body].am.all[x,y,:]

            ax[0].plot(w,am,'x-',label='Component (' + str(x) + ', ' + str(y) + ')')
            ax[1].plot(w,rd,'x-',label='Component (' + str(x) + ', ' + str(y) + ')')
            
            ax[0].legend(loc=0)
                
            plt.show()