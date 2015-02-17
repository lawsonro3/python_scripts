# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:26:30 2013

@author: mlawson
"""
from pylab import *
import flowDistributions as fd

# Determine histrogram, k factor and uBar given in "Wind Data Summary" provided by AWS Truepower, LLC
u = linspace(0,30,100)
weibull = fd.Weibull(k=2.09,uBar=7.3,u=u)
weibull.plotpdf()
weibull.plotcdf()

# Plot the digitized data from "Wind Data Summary" provided by AWS Truepower, LLC
#digitizedData_weib = loadtxt('/home/mlawson/Dropbox/NREL/projects/freshwaterFOA/windData/weibData.dat')
#digitizedData_hist = loadtxt('/home/mlawson/Dropbox/NREL/projects/freshwaterFOA/windData/expData.dat')
#plot(digitizedData_weib[:,0],digitizedData_weib[:,1])
#plot(digitizedData_hist[:,0],digitizedData_hist[:,1])