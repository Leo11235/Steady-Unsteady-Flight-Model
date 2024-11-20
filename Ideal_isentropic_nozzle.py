'''
seeking to model a steady-state isentropic nozzle
Assumptions made:
1. The nozzle flow is one-dimensional; boundary layer flow effects are neglected
2. The gaseous mixture is homogeneous
3. The nozzle gas is in frozen equilibrium, hence its composition does not change across the nozzle axis
4. The fluid obeys the ideal gas law and is calorically perfect, hence the heat capacity ratio is constant
5. The nozzle walls are adiabatic
6. The inlet gas velocity is negligible compared to the exit has velocity

inputs:
Tc: chamber temperature
Wc: chamber gas molar weight
y: heat capacity ratio

we want to calcualte:
At: nozzle throat area
Me: nozzle exit mach number
Te: nozzle gas exit temperature
Ve: nozzle gas exit velocity
F: thrust produced by nozzle
Isp: specific impulse of the engine
'''

class Nozzle():
    def __init__(self, chamber_temperature, molar_weight, heat_capacity_ratio):
        # INPUTS
        # we will use these to calculate the perfect isentropic nozzle
        self.Tc = chamber_temperature
        self.Wc = molar_weight
        self.y = heat_capacity_ratio
