import numpy as np

class SkinFriction(object):
    def model(self, Re) : pass

class LogLaw(SkinFriction):
    def algorithm(self, Re):
        Cf = 0.0576*Re**(-1/5)
        return Cf

class Schlichting(SkinFriction):
    def algorithm(self, Re):
        Cf = (3*np.log10(self.Re) - 0.65)**(-2.3)
        return Cf

class SkinFrictioCalc(object):
    def __init__(self,model=LogLaw())
        self.model = model

    def __init__(self,Re=None,model="power_law"):

        if Re is None:
            print('Re must be provided')
            raise

        else:
            self.Re = Re

        self.model = model

        if self.model=='power_law':

            if self.Re < 5e5 or self.Re > 1e7:
                print('The power law model is only valid for 5e5 < Re <1e7')
                raise

            else:
                self.Cf = 0.0576*Re**(-1/5)

    Schlichting(self):

        if self.Re > 1e9:
            print('The Schlichting model is only valid for < 1e9')
            raise

        else:
            self.Cf = (3*np.log10(self.Re) - 0.65)**(-2.3)
