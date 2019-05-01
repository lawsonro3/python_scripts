'''This object allows skin friction and wall sheer stress to be calculated
'''
import numpy as np

def main():
    '''
    Simple tutorial on how to use this to calculate y+
    '''
    cf = SkinFrictioCalc(Re=1e6, rho=1.2, u=80.0, nu=1.5e-5, y_plus=1.0,model=PowerLaw())

    cf = SkinFrictioCalcMulti(Re=3e6, rho=1.2, u=60.0, nu=1.5e-5, y_plus=1.0)
    cf.calculate(model=PowerLaw())
    cf.calculate(model=PowerLawExp())
    cf.calculate(model=Schlichting())

class SkinFrictionModel(object):
    def __int__(self, Re, rho, u, nu, y_plus):
        self.testVar = 'banana'


    def _cf_calculate(self):
        pass

    def _tau_wall_calculate(self):
        self.tau_wall = 1.0/2.0*self.Cf*self.rho*self.u**2.0
        print('\ttau wall:',self.tau_wall)

    def _u_star_calculate(self):
        self.u_star = np.sqrt(self.tau_wall/self.rho)
        print('\tu*:',self.u_star)

    def _y_calculate(self):
        '''Calculate y for y+ input
        '''
        self.y = self.y_plus*self.nu/self.u_star
        print('\ty @ y+ =',self.y_plus,":",self.y)

class PowerLaw(SkinFrictionModel):
    def _cf_calculate(self):
        self.Cf = 0.0576*self.Re**(-1/5)
        print('1/7th Power Law\n\tCf:',self.Cf)

class PowerLawExp(SkinFrictionModel):
    '''1/7 power law with experimental calibration. Equation 21.12 from
    Schlichting, Hermann (1979), Boundary Layer Theory, ISBN 0-07-055334-3,
    7th Edition.
    '''
    def _cf_calculate(self):
        self.Cf = 0.0592*self.Re**(-1/5)
        print('1/7th Power Law Exp\n\tCf:',self.Cf)

class Schlichting(SkinFrictionModel):
    '''Schlichting, Hermann (1979), Boundary Layer Theory, ISBN 0-07-055334-3,
    7th Edition. Equation 21.16
    '''
    def _cf_calculate(self):
        self.Cf = (2*np.log10(self.Re) - 0.65)**(-2.3)
        print('Schlichting\n\tCf:',self.Cf)

class SkinFrictioCalc(object):
    def __init__(self, Re, rho, u, nu, y_plus, model):
        self.model = model
        self.model.u = u
        self.model.rho = rho
        self.model.nu = nu
        self.model.y_plus = y_plus
        self.model.Re = Re

        self._calculate()

    def _calculate(self):
        self.model._cf_calculate()
        self.model._tau_wall_calculate()
        self.model._u_star_calculate()
        self.model._y_calculate()

    def switch_model(self,model):
        self.model = model
        self.model._calculate()

class SkinFrictioCalcMulti(object):
    def __init__(self, Re, rho, u, nu, y_plus):
            self.u = u
            self.rho = rho
            self.nu = nu
            self.y_plus = y_plus
            self.Re = Re

            self.models = []
            self.model = None

    # def _calculate(self):
    #     for m in self.models:
    #         m._cf_calculate()
    #         m._tau_wall_calculate()
    #         m._u_star_calculate()
    #         m._y_calculate()

    def calculate(self,model):
        self.models.append(model)
        self.models[-1].u = self.u
        self.models[-1].rho = self.rho
        self.models[-1].nu = self.nu
        self.models[-1].y_plus = self.y_plus
        self.models[-1].Re = self.Re

        self.models[-1]._cf_calculate()
        self.models[-1]._tau_wall_calculate()
        self.models[-1]._u_star_calculate()
        self.models[-1]._y_calculate()

# Tutorial
if __name__ == '__main__':
    main()
