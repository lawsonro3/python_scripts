"""
Created on Mon Jul  2 09:32:52 2012

@author: mlawson

This simple program determines the optimal blade shape using Betz theory and blade
element theory. This code is based on section "3.7 - Blade Shape fo Ideal Rotor
without Wake Rotation" in the text book "Wind Energy Explained". Assumptions are:
    1. No wake rotation
    2. No drag
    3. No tip losses
    4. Axial induction factor is 1/3

R       = rotor radius
B       = number of blades
Cl      = lift coeff
tsr     = tip speed ratio
tsr_r   = local tsr at each radial blade station
Cl_Cd_max_aoa = angle of attack at which the maximum lift to drag rario occurs
nElem   = number of elements along the bladel
rStart  = r/R where the blade starts
rEnd    = r/R where the blad
phi     = bla de relative wind or water flow angel
c       = blade chord
c_r     = chord/radius
"""
from pylab import *
interactive(True)

# Define the rotor specs
class Blade(object):
    def __init__(self,R=100.0,B=3.0,Cl=1.0,Cl_Cd_max_aoa=7.0,tsr=7.0,nElem=10,rStart=.1,rEnd=1):
        # User specified or default values        
        self.R = R
        self.B = B
        self.Cl = Cl
        self.tsr = tsr
        self.Cl_Cd_max_aoa = Cl_Cd_max_aoa*pi/180.0 # convert to radians
        
        # Calculated variables
        self.r_R = linspace(rStart,rEnd,nElem)
        self.r = self.r_R*R
        self.tsr_r = self.tsr*self.r_R
        self.phi = array([])      
        for i in range(size(self.tsr_r)):
            self.phi = append(self.phi,math.atan(2.0/3.0/self.tsr_r[i]))
        self.phi_degrees = self.phi*180.0/pi
        self.twist = self.phi - self.Cl_Cd_max_aoa
        self.twist_degrees = self.twist*180.0/pi
        self.c = 8.0*pi*self.r*sin(self.phi)/(3.0*self.B*self.Cl*self.tsr_r)
        self.c_R = self.c/self.R

    def plotBlade(self):
        '''
        	Plots the blade shape
        '''
        figure()
        plot(self.r_R,self.c_R)
        xlabel('r/R')
        ylabel('c/R')
        
        figure()
        plot(self.r_R,self.twist_degrees)
        xlabel('r/R')
        ylabel('Blade twist angle (deg)')
        
    def __str__(self):
        ''' Set out output of the print command'''
        for attr, value in self.__dict__.iteritems():
            print attr, value
        return ''
    def __repr__(self):
        ''' Set out output of the function'''
        for attr, value in self.__dict__.iteritems():
            print attr, value
        return ''
