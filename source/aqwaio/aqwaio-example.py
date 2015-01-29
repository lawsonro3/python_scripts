# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import aqwaio as aio
import matplotlib.pyplot as plt

plt.close('all')

aq = aio.AqwaOutput('/Users/mlawson/Applications/python-scripts/source/aqwaio',outFile='aqwa-data.lis')
aq.plotAddedMassAndDamping(body=0)
aq.plotAddedMassAndDamping(body=1)
aq.plotAddedMassAndDamping(body=2)
aq.writeWecSimHydroData(bodyNumber=1)

