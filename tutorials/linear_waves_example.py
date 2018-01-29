"""
Created on Tue Apr 24 13:49:56 2012

@author: mlawson

Example of how to use the linearWaves Library
"""

from linearWaves import *

# Calculate stuff
wave = CalculateWaveProperities(L=100.,T=8.,H=1.,d=25.)

# Plot stuff
PlotSurfaceElevation(wv=wave)
#PlotParticleDisplacement(wv=wave)
#PlotParticleVelocity(wv=wave)
PlotPressure(wv=wave)
#PlotVelocityPotential(wv=wave)
#PlotParticleAcceleration(wv=wave)