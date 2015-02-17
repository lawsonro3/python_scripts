""""
Created on Wed Apr 25 18:04:08 2012

@author: mlawson
"""
import rainflow
from filters import *
from pylab import *
from scipy.io import loadmat
interactive(True)

dataFile = ''
adcpData = loadmat(dataFile)
V = adcpData['V_nogap']
#samples = linspace(0,size(V,0),size(V,0))
#timeMin = samples*15.
#timeSec = samples*15.*60.
#timeHr = samples/4.Argonne National Laboratory
#timeDay = samples/4./24.type(nonzero(rain.peaksValue==rain.maxPeakVal))

testData = V[0:200,84]
testData,testDataShift = LowPassFilter(testData)
rain = rainflow.CreateData(testData)
rain = rainflow.Count(rain)
rainflow.Hist(rain,10)