''' 
In the original steady model, we assumed a constant mass flow of oxidizer
from the tank to the combustion chamber. 

In this version, mass flow rate depends on the pressure and mass in the tank.
There are several ways to do this, we'll have to decide how to best model this. 
'''

from math import sqrt

class OxidizerTank:
    def __init__(self, tank_volume, tank_pressure, tank_mass, fuel_density, fuel_compressibility, discharge_coefficient, orifice_size):
        self.tank_volume = tank_volume
        self.tank_initial_pressure = tank_pressure
        self.tank_pressure = tank_pressure
        self.tank_mass = tank_mass
        self.fuel_density = fuel_density
        self.fuel_compressibility = fuel_compressibility
        self.discharge_coefficient = discharge_coefficient
        self.orifice_size = orifice_size

    def mrate(self):
        mass_flow_rate = self.discharge_coefficient*sqrt((2*self.tank_pressure)/(self.fuel_density))
        return mass_flow_rate