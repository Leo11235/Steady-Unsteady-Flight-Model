'''
Seeking to model the combustion process in CV2 using the 
semi-empirical hybrid fuel grain regression rate theory. 
Assumptions include:
1. Oxidizer mass flow rate from CV1 to CV2 is constant
2. Combustion chamber pressure is temporally and axially constant
3. Fuel grain regression is axially uniform
4. Combustion gases obey the ideal gas law
5. Combustion chamber walls are adiabatic
'''

from math import pi, sqrt

class CombustionChamber():
    def __init__(self, chamber_pressure, chamber_temperature, chamber_gas_molar_weight, heat_capacity_ratio, regression_rate_scaling_coefficient, regression_rate_exponent, fuel_internal_radius, fuel_external_radius, fuel_length, oxidizer_mass_flow_rate, fuel_grain_density):
        # CONSTANTS
        self.Pc = chamber_pressure
        self.Tc = chamber_temperature
        self.Wc = chamber_gas_molar_weight # NOT THAT # wait maybe that 
        self.y = heat_capacity_ratio
        self.a = regression_rate_scaling_coefficient
        self.n = regression_rate_exponent
        self.Ri = fuel_internal_radius
        self.Re = fuel_external_radius
        self.Lf = fuel_length
        self.MFRo = oxidizer_mass_flow_rate
        self.p = fuel_grain_density

        # MORE CONSTANTS
        self.N = 2 * self.n + 1

        # VARIABLES
        # fuel grain burntime
        self.Tburn = (1 / (self.a * self.N)) * (pi / self.MFRo) ** self.n * (self.Re ** self.N - self.Ri ** 2)
        # fuel grain recession rate (derivative of Ri, separating for variables & integrating in time), not correct, want Ri as a function of time
        self.dRi = (self.a * self.N * (self.MFRo / pi) ** self.n * self.time + rf0 ** self.N) ** (1 / self.N)
        self.dRi = self.a * (self.MFRo / (pi * self.Ri ** 2)) ** self.n
        # remainig fuel grain mass
        self.Mfuel = pi * self.p * self.Lf * (self.Re ** 2 - self.Ri ** 2)
        # fuel mass flow rate, not correct because of dRi
        self.MFRf = 2 * pi * self.p * self.Lf * self.Ri * self.dRi

        # ACTUALLY RELEVANT EQUATIONS
        # fuel grain mass flow rate over the burntime, MFRfavg (eqn 2.6, 2.7)
        # oxidizer-to-fuel ratio, OF (eqn 2.8)
        