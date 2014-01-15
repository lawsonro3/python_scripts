"""
Created on Tue Apr 24 13:49:56 2012

@author: mlawson

Examplt of how to use the linearWaves Library
"""

from linearWaves import *

# Initialize wave properities class
wave = WaveProperties()

# Calculate stuff
wave = CalculateWaveProperities(wv=wave,L=10.,T=2.,H=1.,d=25.)

# Plot stuff
PlotSurfaceElevation(wv=wave)
PlotParticleDisplacement(wv=wave)
PlotParticleVelocity(wv=wave)
PlotPressure(wv=wave)
PlotVelocityPotential(wv=wave)
PlotParticleAcceleration(wv=wave)