'''
These classes define the properities of various materials

@author Michael Lawson
'''
class Fluid(object):
    def __init__(self):
        self.units = {'rho':'kg/m^3',
                     'nu':'m^2/s',
                     'mu':'N-s/m^2'}

class Solid(object):
    def __init__(self):
        self.units = {'rho':'kg/m^3',
                 'E':'Pa',
                 'poissonsRatio':'-',
                 'shearModulus':'Pa',
                 'sigUltimateTensile':'Pa',
                 'sigYield':'Pa'}

class Seawater(Fluid):
    def __init__(self):
        self.rho = 1025.0
        self.nu = 1.0e-6
        self.mu = self.nu*self.rho

class Water(Fluid):
    def __init__(self):    
        self.rho = 1000.0
        self.nu = 1.0e-6
        self.mu = self.nu*self.rho

class Air(Fluid):
    def __init__(self):    
        self.rho = 1.2
        self.nu = 1.0e-5
        self.mu = self.nu*self.rho
    
class Steel(Solid):
    '''
    Steel(grade=Type)
        Available types:
            ASTM-A36
            API 2H Grade 50
    '''
    def __init__(self,grade='ASTM-A36'):
        self.grade = grade        
        if grade is 'API 2H Grade 52': # Need to fill in this option
            self.rho = 'None'
            self.sigYield = 358527378.5
            self.sigUltimateTensile = 455053980.5
            self.E = 'None'
            self.poissonsRatio = 'None'
            self.shearModulus = 'None'
        elif grade is 'API 5L Grade X52': # Need to fill in this option
            self.rho = 'None'
            self.sigYield = 'None'
            self.sigUltimateTensile = 'None'
            self.E = 'None'
            self.poissonsRatio = 'None'
            self.shearModulus = 'None'
        elif grade is 'ASTM-A36':
            self.rho = 7850.0
            self.sigYield = 250.0e6
            self.sigUltimateTensile = 400.036
            self.E = 200.0e9
            self.poissonsRatio = 0.260
            self.shearModulus = 79.3e9
            self.grade = grade

class PolyesterMooringLine(Solid):
    '''
    Approx values for polyester mooring line
    '''
    def __init__(self):
        self.rho = 1380.0
        self.sigYield = 260.0e6
        self.sigUltimateTensile = 426.0e6
        self.E = 14.0e9
    
class Iron(Solid):
    '''
    Nominal values from wikipedia
    '''
    def __init__(self):
        self.rho = 7874.0
        self.E = 211.0e9
        self.sigBreak = 140.0e6

class Concrete(Solid):
    '''
    Nominal values from wikipedia
    '''
    def __init__(self):
        self.rho =  2400.0