# using the various aspects of the rocket previously calculated in calculations.py, we can now obtain all the forces acting on our rocket. 
# we use this to find some information about its ascent by creating a rocket class. 

from variable_initialization import constants_dict
from calculations import calculate_air_density
from math import pi, cos

class Rocket:
    # initialize rocket properties
    def __init__(self, rocket_inputs, rocket_parameters):
        self.num_of_timesteps = constants_dict["number of timesteps"]
        self.timestep_length = rocket_parameters["burntime"] / self.num_of_timesteps

        self.launch_angle = rocket_inputs["launch angle"]
        self.wet_mass = rocket_parameters["wet mass"] # total initial mass

        self.drag_force = 0 # because v=0
        self.grav_force = self.wet_mass * constants_dict["sea level gravity"] * (constants_dict["earth radius"] / (constants_dict["earth radius"] + rocket_inputs["launch site altitude"])) ** 2
        self.thrust_force = rocket_parameters["thrust"]
        self.net_force = self.thrust_force - self.grav_force * cos(self.launch_angle)

        self.vertical_position = rocket_inputs["launch site altitude"]
        self.vertical_velocity = 0
        self.vertial_acceleration = self.net_force / self.wet_mass

    # update values for each timestep, returns a dictionary where keys have list values, each i-th item in a list is that key's value for the i-th timestep
    def calculate_kinematics(self, rocket_inputs, rocket_parameters, constants_dict):
        '''
        The rocket's ascent is split into two parts: during the burn, while it has a force from thrust, and after the burn, where thrust = 0 but v!=0
        We want to calculate apogee, so this sequence ends when v=0
        So this function calls its two helper functions, then unifies their outputs into one
        '''
        burn_dict = during_burn_kinematics(self, rocket_inputs, rocket_parameters, constants_dict)
        ascent_dict = post_burn_kinematics(self, rocket_inputs, rocket_parameters, constants_dict, burn_dict)

        rocket_parameters["reached apogee"] = ascent_dict["position"][-1] # this is how high the rocket reaches

        return burn_dict, ascent_dict



def during_burn_kinematics(self, rocket_inputs, rocket_parameters, constants_dict):
    # initialize unified dictionary with all necessary values for calculating apogee   
    timesteps_dict = {
        "time": [self.timestep_length * i for i in range(self.num_of_timesteps + 1)],
        "drag force": [self.drag_force],
        "gravitational force": [self.grav_force],
        "net force": [self.net_force],
        "mass": [self.wet_mass],
        "acceleration": [self.vertial_acceleration],
        "velocity": [self.vertical_velocity],
        "position": [self.vertical_position],
        }

    # iterate over values to simulate burn
    for i in range(1, self.num_of_timesteps + 1):
        # update forces
        p_air = calculate_air_density(timesteps_dict, constants_dict, i) # calculate air density at current altitude (for drag force)
        new_drag = (1/2) * rocket_inputs["drag coefficient"] * (pi * (rocket_inputs["rocket external diameter"]/2) ** 2) * p_air * timesteps_dict["velocity"][i-1] ** 2
        timesteps_dict["drag force"].append(new_drag)
        new_grav = timesteps_dict["mass"][i-1] * constants_dict["sea level gravity"] * (constants_dict["earth radius"] / (constants_dict["earth radius"] + timesteps_dict["position"][i-1])) ** 2
        timesteps_dict["gravitational force"].append(new_grav)
        new_net_force = self.thrust_force - new_drag - new_grav * cos(self.launch_angle)
        timesteps_dict["net force"].append(new_net_force)
        # update mass
        Mn = rocket_inputs["oxidizer mass flow rate"] + rocket_parameters["average fuel mass flow rate"] # nozzle total propellant mass flow rate
        timesteps_dict["mass"].append(self.wet_mass - Mn * timesteps_dict["time"][i])
        # acceleration
        new_a = timesteps_dict["net force"][i] / timesteps_dict["mass"][i]
        timesteps_dict["acceleration"].append(new_a)
        # velocity
        new_v = timesteps_dict["velocity"][i-1] + timesteps_dict["acceleration"][i] * self.timestep_length
        timesteps_dict["velocity"].append(new_v)
        # position
        new_y = timesteps_dict["position"][i-1] + timesteps_dict["velocity"][i] * self.timestep_length
        timesteps_dict["position"].append(new_y)
    return timesteps_dict



def post_burn_kinematics(self, rocket_inputs, rocket_parameters, constants_dict, burn_dict):
    # set initial values 
    ascent_dict = {
        "time": [burn_dict["time"][-1]], # +self.timestep_length?

        "acceleration": [burn_dict["acceleration"][-1]],
        "velocity": [burn_dict["velocity"][-1]],
        "position": [burn_dict["position"][-1]],

        "mass": burn_dict["mass"][-1], # constant value
        "drag force": [burn_dict["drag force"][-1]],
        "gravitational force": [burn_dict["gravitational force"][-1]],
        "net force": [burn_dict["net force"][-1]],
    }  

    # iterate over values
    i=1
    while ascent_dict["velocity"][-1] > 0: # until reaching apogee,
        # update time
        ascent_dict["time"].append(ascent_dict["time"][-1]+self.timestep_length)
        # update forces
        p_air = calculate_air_density(ascent_dict, constants_dict, i)
        new_drag = (1/2) * rocket_inputs["drag coefficient"] * (pi * (rocket_inputs["rocket external diameter"]/2) ** 2) * p_air * ascent_dict["velocity"][i-1] ** 2
        ascent_dict["drag force"].append(new_drag)
        new_grav = ascent_dict["mass"] * constants_dict["sea level gravity"] * (constants_dict["earth radius"] / (constants_dict["earth radius"] + ascent_dict["position"][i-1])) ** 2
        ascent_dict["gravitational force"].append(new_grav)
        new_net_force = - new_drag - new_grav * cos(rocket_inputs["launch angle"]) # net force is now negative
        ascent_dict["net force"].append(new_net_force)
        # acceleration
        new_a = ascent_dict["net force"][i] / ascent_dict["mass"]
        ascent_dict["acceleration"].append(new_a)
        # velocity
        new_v = ascent_dict["velocity"][i-1] + ascent_dict["acceleration"][i] * self.timestep_length
        ascent_dict["velocity"].append(new_v)
        # position
        new_y = ascent_dict["position"][i-1] + ascent_dict["velocity"][i] * self.timestep_length
        ascent_dict["position"].append(new_y)

        i+=1

    return ascent_dict

# selectively merge the dictionary that keeps information about the burn with the dictionary that keeps info on the ascent
def merge_flight_dictionaries(burn_dict, ascent_dict):

    overall_flight_dict = {
        "time": burn_dict["time"] + ascent_dict["time"],
        "position": burn_dict["position"] + ascent_dict["position"],
        "velocity": burn_dict["velocity"] + ascent_dict["velocity"],
        "acceleration": burn_dict["acceleration"] + ascent_dict["acceleration"],
        "drag force": burn_dict["drag force"] + ascent_dict["drag force"],
        "gravitational force": burn_dict["gravitational force"] + ascent_dict["gravitational force"],
        "net force": burn_dict["net force"] + ascent_dict["net force"],
    }

    return overall_flight_dict
