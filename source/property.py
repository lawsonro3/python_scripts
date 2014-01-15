# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:13:05 2012

A generic properity for a material or a fluid

@author: mlawson
"""
class Property(object):
    value = None
    units = None
    discription = None
    def __init__(self):
        pass
    def __str__(self):
        return 'value = '+str(self.value) \
            + '\nunits = '+str(self.units) \
            + '\ndiscription = '+self.discription               
    def __repr__(self):
        return 'value = '+str(self.value) \
            + '\nunits = '+str(self.units) \
            + '\ndiscription = '+self.discription               