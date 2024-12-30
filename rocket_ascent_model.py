# using the various aspects of the rocket previously calculated in calculations.py, we can now obtain all the forces acting on our rocket. 
# we use this to find some information about its ascent by creating a rocket class. 

from variable_initialization import constants_dict

class Rocket:
    # initialize rocket properties
    def __init__(self, rocket_inputs, rocket_parameters):
        self.num_of_timesteps = constants_dict["number of timesteps"]
        self.timestep_length = rocket_parameters["burntime"] / self.num_of_timesteps
        self.mass = rocket_parameters["wet mass"] # total initial mass
        self.vertical_position = rocket_inputs["launch site altitude"]
        self.vertical_velocity = 0
        self.vertial_acceleration = 0 # GET SUM OF VERTICAL FORCES
        self.drag_force = 0
        self.grav_force = 0
        self.thrust_force = rocket_parameters["thrust"]