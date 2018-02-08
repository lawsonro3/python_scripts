'''
These classes define physical constants. All units are standard metric

@author Michael Lawson
'''
class Constants(object):
    def __init__(self,system=None):
        self.unit_system = 'International M-K-S system'
        self.g = 9.81 # Gravitational acceleration
        self.c_light = 299792458. # Speed of light
        self.c_sount = 340. # Speed of sound
