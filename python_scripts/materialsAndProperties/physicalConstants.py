'''
These classes define physical constants. All units are standard metric

@author Michael Lawson
'''
from property import Property

class Const(object):
    def __init__(self,system=None):
        if system is None:
            self.unitSystem = 'International M-K-S system'
            
            self.g = Property()
            self.g.value = 9.81       
            self.g.units = 'm/s^2'
            self.g.discription = 'Accelleration due to gravity'
    
            self.cLight = Property()
            self.cLight.value = 299792458.
            self.cLight.units = 'm/s'
            self.cLight.discription = 'Speed of light'       
            
            self.cSound = Property()
            self.cSound.value = 340.
            self.cSound.units = 'm/s'
            self.cSound.discription = 'Speed of sound'       
                  
        else:
            self.unitSystem = 'English standard uints'