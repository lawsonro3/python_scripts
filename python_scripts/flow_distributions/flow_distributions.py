# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 14:19:14 2013

@author: mlawson

Wind and water velocity distributions

   
Classses:
    Distributions - Base class that allows flow distribution PDFs and CFDs to be plotted
        plotpdf - function to plot the probability distribution function
            u: bin values
            p: probability
        plodcdf - function to plot the cumulative distribution
            u: bin values
            p: probability
    Raleigh - Rayeigh flow distribution
    Weibull - Weibull flow distribution
    PowerLaw - functionto calculate simple velocity profiles from  power law
"""
from pylab import *
from scipy.special import gamma as gam
interactive(True)

class DistributionPlots(object):
    def plotcdf(self):
        figure('CDF')
        plot(self.u,self.F)
        xlabel('Wind Speed')
        ylabel('Cumulative Distribution Function')
        title(self.distType + ' Distribution')
        
    def plotpdf(self):
        figure('PDF')
        plot(self.u,self.p)
        xlabel('Wind Speed')
        ylabel('Probability Distribution Function')
        title(self.distType + ' Distribution')
               
class Raleigh(DistributionPlots):
    def __init__(self,uBar=None,u=None):
        self.uBar = uBar
        self.distType = 'Raleigh'
        self.u = u
        
        self.p = pi/2.0 * (self.u/self.uBar**2.0) * exp(-pi/4.0*(self.u/self.uBar)**2)
        self.F = 1-exp(-pi/4.0*(self.u/self.uBar)**2)
    
class Weibull(DistributionPlots):
    def __init__(self,k=2.0,uBar=1.0,u=1.0):
        self.distType = 'Weibull'        
        self.k = k
        self.u = u   
        self.uBar = uBar        
        self.c = uBar/math.gamma(1.0+1.0/self.k)
        
        self.p = (self.k/self.c)*(self.u/self.c)**(self.k-1.0)*exp(-(self.u/self.c)**self.k)
        self.F = 1.0-exp(-(self.u/self.c)**self.k)
    
    def update(self,k=None,uBar=None,u=None):
        if k is not None:        
            self.k = k
        else:
            pass
        
        if uBar is not None:
            self.u = u
        else:
            pass
        
        if u is not None:
            self.uBar = uBar
        else:
            pass
        
        self.c = self.uBar/math.gamma(1.0+1.0/self.k)   
         
        self.p = (self.k/self.c)*(self.u/self.c)**(self.k-1.0)*exp(-(self.u/self.c)**self.k)
        self.F = 1.0-exp(-(self.u/self.c)**self.k)
        
class PowerLaw(object):
    def __init__(self,alpha=1.0/7.0,zr=None,ur=None):
        self.alpha = alpha
        self.zr = zr
        self.ur = ur

    def calc_ux(self,zx=None):
        self.zx = zx
        self.ux = self.ur*(self.zx/self.zr)**self.alpha
    
    def calc_zx(self,ux=None):
        self.ux = ux
        self.zx = self.zr*(self.ux/self.ur)**(1/self.alpha)
        
    def plot(self):
        fig_PowerLaw = figure('Power Law')
        plot(ur,zr,'xk',markersize=10,mew=5,label='reference point')
        plot(ux,xz,'xr',markersize=10,mew=5,lavel='reference point')
        legend()