# perform calculations used in the rest of the simulation

from math import pi, sqrt, e

def CV2_calculations(rocket_inputs, rocket_parameters): # control volume 2 = combustion chamber
    rocket_inputs["augmented regression rate exponent"] = calculate_N(rocket_inputs, rocket_parameters)
    rocket_parameters["average oxidizer to fuel ratio"] = calculate_OF(rocket_inputs, rocket_parameters)
    rocket_parameters["fuel mass"] = calculate_fuel_mass(rocket_inputs, rocket_parameters)

# PROPEP calculations fall here

def CV3_calculations(rocket_inputs, rocket_parameters, constants_dict): # CV3 = nozzle
    rocket_parameters["average fuel mass flow rate"] = calculate_Mf(rocket_inputs, rocket_parameters)
    rocket_parameters["total propellant mass flow rate"] = rocket_inputs["oxidizer mass flow rate"] + rocket_parameters["average fuel mass flow rate"]
    rocket_parameters["burntime"] = calculate_Tburn(rocket_inputs, rocket_parameters)
    rocket_parameters["nozzle throat area"] = calculate_At(rocket_inputs, rocket_parameters, constants_dict)
    rocket_parameters["nozzle throat radius"] = calculate_Rt(rocket_inputs, rocket_parameters)

    #rocket_parameters["nozzle gas exit pressure"] = calculate_Pe(rocket_inputs, rocket_parameters, constants_dict)
    # playing with exit pressure gives better results, very strange
    # in an ideal nozzle, exit pressure should be the same as 1 atm
    rocket_parameters["nozzle gas exit pressure"] = 101325 * 0.95926

    rocket_parameters["nozzle gas exit mach number"] = calculate_Me(rocket_inputs, rocket_parameters)
    rocket_parameters["nozzle exit area"] = calculate_Ae(rocket_inputs, rocket_parameters)
    rocket_parameters["nozzle exit radius"] = calculate_Re(rocket_inputs, rocket_parameters)
    rocket_parameters["nozzle gas exit temperature"] = calculate_Te(rocket_inputs, rocket_parameters)
    rocket_parameters["nozzle gas exit velocity"] = calculate_Ve(rocket_inputs, rocket_parameters, constants_dict)
    rocket_parameters["thrust"] = calculate_F(rocket_inputs, rocket_parameters, constants_dict)
    rocket_parameters["Isp"] = calculate_Isp(rocket_inputs, rocket_parameters, constants_dict)
    rocket_parameters["total impulse"] = calculate_Ns(rocket_inputs, rocket_parameters)
    rocket_parameters["wet mass"] = calculate_Mw(rocket_inputs, rocket_parameters)
    rocket_parameters["thrust to weight ratio"] = calculate_TtW(rocket_inputs, rocket_parameters, constants_dict)


### CALCULATIONS BELOW ###

# augmented regression rate exponent
def calculate_N(rocket_inputs, rocket_parameters):
    n = rocket_inputs["regression rate exponent"]
    
    N = 2 * n + 1
    return N

# average oxidizer to fuel ratio
def calculate_OF(rocket_inputs, rocket_parameters):
    Mo = rocket_inputs["oxidizer mass flow rate"]
    a = rocket_inputs["regression rate scaling coefficient"]
    n = rocket_inputs["regression rate exponent"]
    N = rocket_inputs["augmented regression rate exponent"]
    p = rocket_inputs["fuel grain density"]
    Lf = rocket_inputs["fuel length"]
    Re = rocket_inputs["fuel external radius"]
    Ri0 = rocket_parameters["initial internal fuel radius"]
    
    OF = (1 / (a * N * p * Lf)) * (Mo / pi) ** (1 - n) * ((Re ** N - Ri0 ** N) / (Re ** 2 - Ri0 ** 2))
    return OF 

def calculate_fuel_mass(rocket_inputs, rocket_parameters):
    Lf = rocket_inputs["fuel length"]
    Ri0 = rocket_parameters["initial internal fuel radius"]
    Re = rocket_inputs["fuel external radius"] 
    p = rocket_inputs["fuel grain density"]

    return pi * Lf * (Re ** 2 - Ri0 ** 2) * p

# average mass flow rate
def calculate_Mf(rocket_inputs, rocket_parameters):
    Mo = rocket_inputs["oxidizer mass flow rate"]
    a = rocket_inputs["regression rate scaling coefficient"]
    n = rocket_inputs["regression rate exponent"]
    N = rocket_inputs["augmented regression rate exponent"]
    p = rocket_inputs["fuel grain density"]
    Lf = rocket_inputs["fuel length"]
    Re = rocket_inputs["fuel external radius"] 
    Ri0 = rocket_parameters["initial internal fuel radius"]

    Mf = a * pi * N * p * Lf * ((Re ** 2 - Ri0 ** 2) / (Re ** N - Ri0 ** N)) * (Mo / pi) ** n
    return Mf

# burntime
def calculate_Tburn(rocket_inputs, rocket_parameters):
    Mo = rocket_inputs["oxidizer mass flow rate"]
    a = rocket_inputs["regression rate scaling coefficient"]
    n = rocket_inputs["regression rate exponent"]
    N = rocket_inputs["augmented regression rate exponent"]
    Re = rocket_inputs["fuel external radius"]
    Ri0 = rocket_parameters["initial internal fuel radius"]

    Tburn = (1 / (a * N * (Mo / pi) ** n)) * (Re ** N - Ri0 ** N)
    return Tburn

# calculate ideal nozzle throat area
def calculate_At(rocket_inputs, rocket_parameters, constants_dict):
    Mn = rocket_inputs["oxidizer mass flow rate"] + rocket_parameters["average fuel mass flow rate"] # nozzle total propellant mass flow rate
    Pc = rocket_inputs["chamber pressure"] 
    Mt = 1 # nozzle throat mach number, this is always the case for a converging diverging engine
    Ru = constants_dict["universal gas constant"]
    Tc = rocket_parameters["chamber temperature"]
    gamma = rocket_parameters["heat capacity ratio"]
    Wc = rocket_parameters["chamber gas molar weight"]

    At = (Mn / (Pc * Mt)) * sqrt((Ru * Tc) / (gamma * Wc)) * (1 + ((gamma - 1) / 2) * Mt ** 2) ** ((gamma + 1) / (2 * (gamma - 1)))
    return At

def calculate_Rt(rocket_inputs, rocket_parameters):
    At = rocket_parameters["nozzle throat area"]

    return sqrt(At / pi)

# calculate nozzle gas exit pressure
# this one is a bit unique, no actual calculation because we assume that the nozzle is perfectly expanded
# this means the exit pressure matches ambient pressure
def calculate_Pe(rocket_inputs, rocket_parameters, constants_dict):
    Pinf = constants_dict["ambient sea level atmospheric pressure"]
    #Pe = Pinf * .53102492 # no clue why but multiplying by this number solves all my problems
    Pe = Pinf
    return Pe

# nozzle gas exit mach number
def calculate_Me(rocket_inputs, rocket_parameters):
    gamma = rocket_parameters["heat capacity ratio"]
    Pc = rocket_inputs["chamber pressure"]
    Pe = rocket_parameters["nozzle gas exit pressure"]

    Me = sqrt((2 / (gamma - 1)) * ((Pc / (Pe)) ** ((gamma - 1) / (gamma)) - 1)) 
    return Me

# ideal nozzle exit area
def calculate_Ae(rocket_inputs, rocket_parameters):
    At = rocket_parameters["nozzle throat area"]
    Me = rocket_parameters["nozzle gas exit mach number"]
    gamma = rocket_parameters["heat capacity ratio"]

    Ae = (At / Me) * ((2/(gamma+1))*(1 + ((gamma-1)/(2)) * Me ** 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
    return Ae

# nozzle exit radius
def calculate_Re(rocket_inputs, rocket_parameters):
    Ae = rocket_parameters["nozzle exit area"]

    return sqrt(Ae / pi)

# nozzle gas exit temperature
def calculate_Te(rocket_inputs, rocket_parameters):
    Tc = rocket_parameters["chamber temperature"]
    gamma = rocket_parameters["heat capacity ratio"]
    Me = rocket_parameters["nozzle gas exit mach number"]

    Te = Tc / (1 + ((gamma - 1) / (2)) * Me ** 2)
    return Te

# calculate nozzle gas exit velocity
def calculate_Ve(rocket_inputs, rocket_parameters, constants_dict):
    Me = rocket_parameters["nozzle gas exit mach number"]
    gamma = rocket_parameters["heat capacity ratio"]
    Te = rocket_parameters["nozzle gas exit temperature"]
    Ru = constants_dict["universal gas constant"]
    Wc = rocket_parameters["chamber gas molar weight"]

    Ve = Me * sqrt((gamma * Te * Ru) / Wc)
    return Ve

# calculate thrust
def calculate_F(rocket_inputs, rocket_parameters, constants_dict):
    Mn = rocket_inputs["oxidizer mass flow rate"] + rocket_parameters["average fuel mass flow rate"] # nozzle total propellant mass flow rate
    Ve = rocket_parameters["nozzle gas exit velocity"]
    Pe = rocket_parameters["nozzle gas exit pressure"]
    Pinf = constants_dict["ambient sea level atmospheric pressure"]
    Ae = rocket_parameters["nozzle exit area"]

    F = Mn * Ve + (Pe - Pinf) * Ae
    return F

# engine Isp
def calculate_Isp(rocket_inputs, rocket_parameters, constants_dict):
    F = rocket_parameters["thrust"]
    Mn = rocket_inputs["oxidizer mass flow rate"] + rocket_parameters["average fuel mass flow rate"] # nozzle total propellant mass flow rate
    Gsl = constants_dict["sea level gravity"]

    Isp = F / (Mn * Gsl)
    return Isp

# total impusle
def calculate_Ns(rocket_inputs, rocket_parameters):
    F = rocket_parameters["thrust"]
    t = rocket_parameters["burntime"]

    return F * t

# wet mass
def calculate_Mw(rocket_inputs, rocket_parameters):
    Md = rocket_inputs["dry mass"]
    Mf = rocket_parameters["fuel mass"]
    t = rocket_parameters["burntime"]
    Mo = rocket_inputs["oxidizer mass flow rate"]

    Mw = Md + Mf + Mo * t
    return Mw

# thrust to weight ratio
def calculate_TtW(rocket_inputs, rocket_parameters, constants_dict):
    thrust = rocket_parameters["thrust"]
    mass = rocket_parameters["wet mass"]
    Gsl = constants_dict["sea level gravity"]

    return thrust / (mass * Gsl)

# calculates air density for a given altitude (for rocket_ascent_model.py)
def calculate_air_density(constants_dict, height):
    T0 = constants_dict["sea level temperature"]
    T11 = constants_dict["stratosphere temperature"]
    L = constants_dict["temperature lapse rate in the troposphere"]
    p0 = constants_dict["sea level air density"]
    g = constants_dict["sea level gravity"]
    M = constants_dict["air molar mass"]
    Ru = constants_dict["universal gas constant"]
    if height <= 11000: # if the rocket is in the troposphere
        p_air = p0 * ((T0 - L * height) / T0) ** (g*M/(Ru*L) - 1)
    else: # the rocket is in the stratosphere
        p11 = p0 * (T11 / T0) ** (g*M/(Ru*L) - 1) # pressure at 11000 m
        p_air = p11 * e ** ((11000 - height)*(g*M/(Ru*T11)))

    return p_air

# calculate gravity as a function of height
def calculate_gravity(constants_dict, height):
    R = constants_dict["earth radius"]
    return constants_dict["sea level gravity"] * (R+height)/R