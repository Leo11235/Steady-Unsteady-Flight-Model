from math import e


# calculates air density for a given altitude
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

# returns ambient pressure at rocket_height
def calculate_ambient_pressure(constants_dict, rocket_height):
    return