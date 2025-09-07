# using the various aspects of the rocket previously calculated in calculations.py, we can now obtain all the forces acting on our rocket. 
# we use this to find some information about its ascent by creating a rocket class. 

from steady_backend.calculations import calculate_air_density, calculate_gravity
from math import pi, cos, radians

# takes as input details about the rocket, outputs a dict with flight information
def simulate_rocket_ascent(rocket_inputs, rocket_parameters, constants_dict):
    helper_dict = { # time & constants
        "drag constant" : 0.5 * rocket_inputs["drag coefficient"] * pi * (rocket_inputs["rocket external diameter"]/2) ** 2,
        "cos angle" : cos(radians(rocket_inputs["launch angle"])),
        "timestep length": rocket_parameters["burntime"] / constants_dict["number of timesteps in burn"], 
        "vertical thrust": rocket_parameters["thrust"] * cos(radians(rocket_inputs["launch angle"]))
    }
    # calculate ascent
    burn_dict = calculate_burn_kinematics(rocket_inputs, rocket_parameters, constants_dict, helper_dict)
    ascent_dict = calculate_post_burn_kinematics(rocket_inputs, rocket_parameters, constants_dict, burn_dict, helper_dict)
    flight_dict = merge_flight_dictionaries(burn_dict, ascent_dict)
    return flight_dict
    
    
# calculate burn kinematics
def calculate_burn_kinematics(rocket_inputs, rocket_parameters, constants_dict, helper_dict):
    dt = helper_dict["timestep length"]
    init_grav_force = rocket_parameters["wet mass"] * calculate_gravity(constants_dict, rocket_inputs["launch site altitude"]) * (-1) # *-1 because gravity points down
    init_net_force = rocket_parameters["thrust"] * helper_dict["cos angle"] + init_grav_force
    burn_dict = { # create flight dict, set initial parameters (time, thrust, and mass remain the same or change predictably and we can set them from the getgo)
        "time": [helper_dict["timestep length"] * i for i in range(constants_dict["number of timesteps in burn"]+1)], # used to keep a record of exactly when everything happens
        "thrust": [helper_dict["vertical thrust"] for i in range(constants_dict["number of timesteps in burn"]+1)], # thrust remains constant throughout burn
        "drag force": [0],
        "grav force": [init_grav_force],
        "net force": [init_net_force],
        "mass": [rocket_parameters["wet mass"] - i * rocket_parameters["total propellant mass flow rate"] * dt for i in range(constants_dict["number of timesteps in burn"]+1)],
        "acceleration": [init_net_force / rocket_parameters["wet mass"]],
        "velocity": [0],
        "altitude": [rocket_inputs["launch site altitude"]],
    }

    # iterate over values to simulate burn
    for i in range(1, constants_dict["number of timesteps in burn"]+1):
        # set variables
        v = burn_dict["velocity"][-1]
        h = burn_dict["altitude"][-1]
        m = burn_dict["mass"][i]
        # update vertical forces
        thrust = helper_dict["vertical thrust"]
        drag = helper_dict["drag constant"] * calculate_air_density(rocket_inputs, constants_dict, h) * abs(v) * (-1) * v # no need to multiply by cos(angle) since v is already the vertical value
        grav = (-1) * m * calculate_gravity(constants_dict, h)
        net_force = thrust + drag + grav
        # calculate kinematics
        new_a = net_force/burn_dict["mass"][i]
        new_v = burn_dict["velocity"][-1] + new_a * dt
        new_h = burn_dict["altitude"][-1] + new_v * dt
        # append new values
        burn_dict["drag force"].append(drag)
        burn_dict["grav force"].append(grav)
        burn_dict["net force"].append(net_force)
        burn_dict["acceleration"].append(new_a)
        burn_dict["velocity"].append(new_v)
        burn_dict["altitude"].append(new_h)
    return burn_dict
    
    
# iterate over values to simulate postburn kinematics
def calculate_post_burn_kinematics(rocket_inputs, rocket_parameters, constants_dict, burn_dict, helper_dict):
    dt = helper_dict["timestep length"]
    ascent_dict = { # create flight dict, set initial parameters (time, thrust, and mass remain the same or change predictably and we can set them from the getgo)
        # first set of values is 'stolen' from burn_dict, will get deleted later
        "time": [burn_dict["time"][-1]], # used to keep a record of exactly when everything happens
        "thrust": [0], # thrust remains constant throughout burn
        "drag force": [burn_dict["drag force"][-1]],
        "grav force": [burn_dict["grav force"][-1]],
        "net force": [burn_dict["net force"][-1]],
        "mass": [burn_dict["mass"][-1]], # stays constant throughout this phase
        "acceleration": [burn_dict["acceleration"][-1]],
        "velocity": [burn_dict["velocity"][-1]],
        "altitude": [burn_dict["altitude"][-1]],
    }
    
    m = ascent_dict["mass"][0]
    while ascent_dict["velocity"][-1] > 0: # until reaching apogee, 
        # update time, append constants
        ascent_dict["time"].append(ascent_dict["time"][-1]+dt)
        ascent_dict["mass"].append(ascent_dict["mass"][0])
        ascent_dict["thrust"].append(0)
        # set variables
        v = ascent_dict["velocity"][-1]
        h = ascent_dict["altitude"][-1]
        # update vertical forces
        drag = helper_dict["drag constant"] * calculate_air_density(rocket_inputs, constants_dict, h) * abs(v) * (-1) * v # no need to multiply by cos(angle) since v is already the vertical value
        grav = (-1) * m * calculate_gravity(constants_dict, h)
        net_force = drag + grav
        # calculate kinematics
        new_a = net_force/m
        new_v = ascent_dict["velocity"][-1] + new_a * dt
        new_h = ascent_dict["altitude"][-1] + new_v * dt
        # append new values
        ascent_dict["drag force"].append(drag)
        ascent_dict["grav force"].append(grav)
        ascent_dict["net force"].append(net_force)
        ascent_dict["acceleration"].append(new_a)
        ascent_dict["velocity"].append(new_v)
        ascent_dict["altitude"].append(new_h)
    return ascent_dict
    
# take burn & ascent dicts, merge them while deleting the first column of information from ascent_dict
def merge_flight_dictionaries(burn_dict, ascent_dict):
    flight_dict = {}
    for key in burn_dict:
        flight_dict[key] = burn_dict[key] + ascent_dict[key][1:]
    return flight_dict