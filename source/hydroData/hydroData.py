# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:50:35 2015

@author: mlawson
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
plt.interactive(True)

class HydrodynamicCoefficients(object):
    def __init__(self):
        self.all            = np.array([])
        self.inf        = np.array([])
        self.zero       = np.array([])
    
class HydrodynamicExcitation(object):
    def __init__(self):
        self.re             = np.array([])
        self.im             = np.array([])
        self.mag            = np.array([])
        self.phase          = np.array([])
        
class HydrodynamicData(object):
    def __init__(self):
        self.rho            = 1000
        self.g              = 9.81      
        self.files          = {}
        self.nBodies        = 0                             # Number of bodies in simulation
        self.cg             = {}                            # Center of gravity
        self.cb             = {}                            # Center of buoyancy
        self.volDisp        = {}                            # Volume displacement
        self.T              = {}                            # Wave period
        self.w              = {}                            # Wave freq
        self.am             = HydrodynamicCoefficients()    # Added mass 
        self.rd             = HydrodynamicCoefficients()    # Radiation damping
        self.wpArea         = {}                            # Water plane area          
        self.buoyForce      = {}                            # Buoyanch force at equelibrium
        self.k              = {}                            # Hydrostatic stifness matrix
        self.ex             = HydrodynamicExcitation()      # Excitation coeffs
        self.waterDepth     = None                          # Water depth
        self.waveDir        = 0                             # Wave direction
        self.name           = None                          # Name of the body
        
def writeHdf5(data,outFile):

        try:

            import h5py

        except:

            raise Exception('The h5py module must be installed to used the writeHdf5 functionality.')
            
            
            
        with h5py.File(outFile, "w") as f:        

            for i in range(np.size(data.keys())):

                T = f.create_dataset('body' + str(i) + '/sim/T',data=data[i].T)
                T.attrs['units'] = 's'
                T.attrs['description'] = 'Wave periods'
                
                w = f.create_dataset('body' + str(i) + '/sim/w',data=data[i].w)
                w.attrs['units'] = 'rad/s'                
                w.attrs['description'] = 'Wave frequencies'                
                
                k = f.create_dataset('body' + str(i) + '/body/k',data=data[i].k)
                k.attrs['units'] = ''
                k.attrs['description'] = 'Hydrostatic stiffness matrix'  
                
                
                cg = f.create_dataset('body' + str(i) + '/body/cg',data=data[i].cg)
                cg.attrs['units'] = 'm'
                cg.attrs['description'] = 'Center of gravity'  

                cb = f.create_dataset('body' + str(i) + '/body/cb',data=data[i].cb)
                cb.attrs['units'] = 'm'
                cb.attrs['description'] = 'Center of buoyancy'  
                
                exMag = f.create_dataset('body' + str(i) + '/hydro/ex/mag',data=data[i].ex.mag)
                exMag.attrs['units'] = ''
                exMag.attrs['description'] = 'Magnitude of excitation force'  
                
                exPhase = f.create_dataset('body' + str(i) + '/hydro/ex/phase',data=data[i].ex.phase)
                exPhase.attrs['units'] = 'rad'
                exPhase.attrs['description'] = 'Phase angle of exctiation force'  
                
                exRe = f.create_dataset('body' + str(i) + '/hydro/ex/re',data=data[i].ex.re)
                exRe.attrs['units'] = ''
                exRe.attrs['description'] = 'Real component of excitation force'  

                exIm = f.create_dataset('body' + str(i) + '/hydro/ex/im',data=data[i].ex.im)
                exIm.attrs['units'] = ''
                exIm.attrs['description'] = 'Imaginary component of excitation force'  

                # Write added mass information                
                amInf = f.create_dataset('body' + str(i) + '/hydro/am/inf',data=data[i].am.infFreq)
                amInf.attrs['units for translational degrees of freedom'] = 'kg'
                amInf.attrs['description'] = 'Infinite frequency added mass'
                
                am = f.create_dataset('body' + str(i) + '/hydro/am/all',data=data[i].am.all)
                am.attrs['units for translational degrees of freedom'] = 'kg'                
                am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
                am.attrs['description'] = 'Added mass. Frequency is the thrid dimension of the data structure.'
                
                for m in xrange(np.shape(data[i].am.all)[0]):
                
                    for n in xrange(np.shape(data[i].am.all)[1]):

                        amComp = f.create_dataset('body' + str(i) + '/hydro/am/comps/' + str(m) + ',' + str(n),data=data[i].am.all[m,n,:])
                        amComp.attrs['units'] = ''
                        amComp.attrs['description'] = 'Added mass components as a function of frequency'

                        radComp = f.create_dataset('body' + str(i) + '/hydro/rd/comps/' + str(m) + ',' + str(n),data=data[i].rd.all[m,n,:])
                        radComp.attrs['units'] = ''
                        radComp.attrs['description'] = 'Radiation damping components as a function of frequency'
                
                rad = f.create_dataset('body' + str(i) + '/hydro/rd/all',data=data[i].rd.all)
                rad.attrs['units'] = ''
                rad.attrs['description'] = 'Radiation damping. Frequency is the thrid dimension of the data structure.'
                
                wDepth = f.create_dataset('body' + str(i) + '/sim/wDepth',data=data[i].waterDepth)
                wDepth.attrs['units'] = 'm'
                wDepth.attrs['description'] = 'Water depth'

                waveHead = f.create_dataset('body' + str(i) + '/sim/wDir',data=data[i].waveDir)
                waveHead.attrs['units'] = 'rad'
                waveHead.attrs['description'] = 'Wave direction'
                
                vol = f.create_dataset('body' + str(i) + '/body/dispVol',data=data[i].volDisp)
                vol.attrs['units'] = 'm^3'
                vol.attrs['description'] = 'Displaced volume'

                
                g = f.create_dataset('body' + str(i) + '/sim/g',data=data[i].g)
                g.attrs['units'] = 'm/s^2'
                g.attrs['description'] = 'Gravitational acceleration'
                
                rho = f.create_dataset('body' + str(i) + '/sim/rho',data=data[i].rho)
                rho.attrs['units'] = 'kg/m^3'
                rho.attrs['description'] = 'Water density'
                
            print 'Wrote HDF5 data to ' + outFile


    
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
            out['fDamping'] = curData.rd.all[:,:,::-1]
            out['fExtRe'] = curData.ex.re[::-1,:].transpose()
            out['fExtIm'] = curData.ex.im[::-1,:].transpose()
            out['fExtMag'] = curData.ex.mag[::-1,:].transpose()
            out['fExtPhase'] = curData.ex.phase[::-1,:].transpose()
            
            outFileName = outFile+'-body' + str(i) +'.mat'
            sio.savemat(outFileName,out)

            print 'Wrote MATLAB output for WEC-Sim to ' + outFileName
            
            
            
def plotAddedMassAndDamping(data,fName,components):
    '''
    Function to plot the diagional component of added mass and raditation
    damping - this could be singificantly improved
    Inputs:
        A list of components to plot
    Outputs:
        None
    '''
    for body in xrange(np.size(data.keys())):
        
        fNameTemp = fName + '-addedMassRadDamping-body' + str(body) + '.ps'
        
        f, ax = plt.subplots(2, sharex=True, figsize=(7,10))
        ax[0].plot()
        ax[0].set_title('Hydrodynamic coefficients for body ' + str(body) + ':' + str(data[body].name))    
        ax[0].set_ylabel('Added mass')
        
        ax[1].plot()

        ax[1].set_xlabel('Wave frequency (rad/s)')
        ax[1].set_ylabel('Radiation damping')
        
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
        f.savefig(fNameTemp,dpi=300)

def plotExcitation(data,fName,components):

    for body in xrange(np.size(data.keys())):
        
        fNameTemp = fName + '-excitation-body' + str(body) + '.ps'
        
        f, ax = plt.subplots(4, sharex=True,figsize=(7,10))
        ax[0].plot()
        ax[0].set_title('Excitation force for body ' + str(body) + ':' + str(data[body].name))    
        ax[0].set_ylabel('Ex force - real')
        
        ax[1].plot()
        ax[1].set_ylabel('Ex force - imaginary')
        
        ax[2].plot()
        ax[2].set_ylabel('Ex force - mag')
        
        ax[3].plot()        
        ax[3].set_xlabel('Wave frequency (rad/s)')        
        ax[3].set_ylabel('Ex force - phase')
        
        for i,comp in enumerate(components):
            
            m = comp
            w = data[body].w
            re = data[body].ex.re[:,m]
            im = data[body].ex.im[:,m]
            mag = data[body].ex.mag[:,m]
            phase = data[body].ex.phase[:,m]

            ax[0].plot(w,re,'x-',label='Component (' + str(m) + ')')
            ax[1].plot(w,im,'x-',label='Component (' + str(m) + ')')
            ax[2].plot(w,mag,'x-',label='Component (' + str(m) + ')')
            ax[3].plot(w,phase,'x-',label='Component (' + str(m) + ')')
            
            ax[0].legend(loc=0)
                
        plt.show()
        f.savefig(fNameTemp,dpi=300)