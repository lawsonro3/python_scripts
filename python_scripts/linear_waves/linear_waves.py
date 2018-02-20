"""
Created on Tue Apr 24 09:29:26 2012

@author: mlawson

Calculate properities of linear waves. Equations from Sarpkaya "Wave Forces on
Offshore Structures", page 113

Variable Definitions
d                           water depth
H                           wave height
L                           wave length
T                           wave period
x=0                         x coord - parallel to the free surface
z=0                         z coord - normal to the free surface
k                           wave number
c                           wave velocity
t                           physical time
rho=1025                    fluid density
cg                          group velocity
P                           energy flux
phi                         velocity potential

How to use
1. Load the
2. Execute the SetWaveProperities to define the wave
3. Then call any of the other functions to calculate data or make plots

"""
from pylab import *
interactive(True)

class WaveProperties(object):
    T = float()
    L = float()
    H = float()
    d = float()
    x = array([])
    z = array([])
    k = float()
    omega = float()
    c = float()
    s = float()
    E = float()
    cg = float()
    P = float()
    rho = float()
    g = float()
    t = array([])
    theta = array([])
    xDisp = array([])
    zDisp = array([])
    u = array([])
    w = array([])
    du_dt = array([])
    dw_dt = array([])

def CalculateWaveProperities(L,T,H,d,tLen=25,xLen=100,zLen=75,rho=1025.0,g=9.81):

    # Load class
    wv = WaveProperties()

    # Period
    wv.T = float(T)

    # Wave length
    wv.L = float(L)

    # Wave height
    wv.H = float(H)

    # Water Depth
    wv.d = float(d)

    # Fluid density
    wv.rho = float(rho)

    # Gravity
    wv.g = float(g)

    # Spatial vectors - this is a shitty way to do a "3D meshgrid"
    xr = linspace(0,L,xLen)
    zr = linspace(-d,0,zLen)
    xr,zr = meshgrid(xr,zr)
    wv.x = repeat(xr,tLen).reshape(zLen,xLen,tLen)
    wv.z = repeat(zr,tLen).reshape(zLen,xLen,tLen)

    # Temporal vector - this is a shitty way to do a "3D meshgrid"
    t = linspace(0,T,tLen)
    wv.t = ones(shape(wv.x))
    for i in range(tLen):
        wv.t[:,:,i] = wv.t[:,:,i]*t[i]

    # Wave number
    wv.k = 2.*pi/wv.L

    # Angular frequency
    wv.omega = 2.*pi/wv.T

    # Wave velocity
    wv.c = wv.omega/wv.k

    # Something ?????
    wv.s = wv.z+wv.d

    # Something ?????
    wv.theta = wv.k*(wv.x-wv.c*wv.t)

    # Average wave energy density
    wv.E = 1./8.*wv.rho*wv.g*wv.H**2.

    # Group velocity
    wv.cg = 1./2.*(1.+2.*wv.k*wv.d/sinh(2*wv.k*wv.d))*wv.c

    # Energy flux
    wv.P = wv.E*wv.cg

    # Velocity potential
    wv.phi = pi*wv.H/wv.k/wv.T*cosh(wv.k*wv.s)/sinh(wv.k*wv.d)*sin(wv.theta)

    # Particle displacement
    wv.xDisp = -wv.H/2.*cosh(wv.k*wv.s)/sinh(wv.k*wv.d)*sin(wv.theta)
    wv.zDisp = wv.H/2.*sinh(wv.k*wv.s)/sinh(wv.k*wv.d)*cos(wv.theta)

    # Surface elevation
    wv.eta = wv.H/2.*cos(wv.theta)

    # Particle velocity
    wv.u = pi*wv.H/wv.T*cosh(wv.k*wv.s)/sinh(wv.k*wv.d)*cos(wv.theta)
    wv.w = pi*wv.H/wv.T*sinh(wv.k*wv.s)/sinh(wv.k*wv.d)*sin(wv.theta)

    # Particle accelleration
    wv.du_dt = 2.*pi**2.*wv.H/wv.T**2.*cosh(wv.k*wv.s)/sinh(wv.k*wv.d)*sin(wv.theta)
    wv.dw_dt = -2.*pi**2.*wv.H/wv.T**2.*sinh(wv.k*wv.s)/sinh(wv.k*wv.d)*cos(wv.theta)

    # Pressure
#    wv.p = -wv.rho*wv.g*wv.z + 1./2.*wv.rho*wv.g*wv.H*cosh(wv.k*wv.s)/cosh(wv.k*wv.d)*cos(wv.theta)
    wv.p = 1./2.*wv.rho*wv.g*wv.H*cosh(wv.k*wv.s)/cosh(wv.k*wv.d)*cos(wv.theta)

    return wv

def PlotSurfaceElevation(wv,xPos=0,zPos=0):
        figure('Surface elevation')
        plot(wv.t[zPos,xPos,:],wv.eta[zPos,xPos,:])
        xlabel('Time (s)')
        ylabel('Surface elevation')

def PlotParticleDisplacement(wv,timeInd=0):
        figure('Particle displacements')
        subplot(2,1,1)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.xDisp[:,:,timeInd],100)
        ylabel('z')
        colorbar()
        title('x particle displacement')
        subplot(2,1,2)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.zDisp[:,:,timeInd],100)
        xlabel('x')
        ylabel('z')
        title('z particle displacement')
        colorbar()

def PlotParticleVelocity(wv,timeInd=0):
        figure('Particle velocity')
        subplot(2,1,1)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.u[:,:,timeInd],100)
        ylabel('z')
        title('u velocity')
        colorbar()
        subplot(2,1,2)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.w[:,:,timeInd],100)
        title('u velocity')
        ylabel('z')
        xlabel('x')
        colorbar()

def PlotParticleAcceleration(wv,timeInd=0):
        figure('Particle acceleration')
        subplot(2,1,1)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.du_dt[:,:,timeInd],100)
        ylabel('z')
        title('u accelleration')
        colorbar()
        subplot(2,1,2)
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.dw_dt[:,:,timeInd],100)
        title('w acceleration')
        ylabel('z')
        xlabel('x')
        colorbar()

def PlotPressure(wv,timeInd=0):
        figure('Pressure')
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.p[:,:,timeInd],100)
        xlabel('x')
        ylabel('z')
        colorbar()

def PlotVelocityPotential(wv,timeInd=0):
        figure('Velocity potential')
        contourf(wv.x[:,:,timeInd],wv.z[:,:,timeInd],wv.p[:,:,timeInd],100)
        xlabel('x')
        ylabel('z')
        colorbar()


if __name__ == '__main__':
    # Calculate stuff
    wave = CalculateWaveProperities(L=100.,T=8.,H=1.,d=25.)

    # Plot stuff
    PlotSurfaceElevation(wv=wave)
    #PlotParticleDisplacement(wv=wave)
    #PlotParticleVelocity(wv=wave)
    PlotPressure(wv=wave)
    #PlotVelocityPotential(wv=wave)
    #PlotParticleAcceleration(wv=wave)
