'''
this class handless the input file and format is in a way
that is nice to use for the rest of the program
'''

import json

# list of required inputs
required_keys = [
    "rocket name",

    # combustion stuff
    "oxidizer mass flow rate",
    "chamber pressure",
    "fuel external diameter",
    "fuel length",

    # rocket info
    "target apogee",
    "launch site altitude",
    "dry mass",
    "rocket external diameter",
    "drag coefficient",
    "launch angle",

    # idk yet
    "prop dry mass",
    "prop diameter",
    "prop length",

    # for PROPEP
    "liquid oxidizer type",
    "oxidizer molecule code",
    "solid fuel type",
    "fuel molecule code",

    # from research
    "fuel grain density",
    "regression rate scaling coefficient",
    "regression rate exponent",
]

# dictionary containing universal constants
# may end up adding lots of random stuff here
constants_dict = {
    # natural constants
    "sea level gravity": 9.80665, # m/s^2
    "universal gas constant": 8.31446261815324, # J / (mol * K)
    "ambient sea level atmospheric pressure": 101325,  # Pa 
    "earth radius": 6378000, # m
    "sea level air density": 1.225, # kg/m^3
    "sea level temperature": 288.15, # K
    "stratosphere temperature": 216.65, # K
    "air molar mass": 0.0289644, # kg/mol
    "temperature lapse rate in the troposphere": 0.0065, # K/m

    # for computation purposes
    "number of timesteps": 100, # number of timesteps in burntime (only accounts for the part of the burn where there is a thrust). Ascent timesteps typically 300-400
    "tolerated apogee difference": 1, # m, allowed space between target and real apogee
}

# main function of this class, takes input file, outputs a dictionary rocket_inputs
# this entire function may get deleted if I can get the frontend to directly send a dictionary to this code
def read_rocket_attributes(json_file_path):
    rocket_inputs = {}

    # read JSON file and populate rocket_inputs
    try:
        with open(json_file_path, 'r') as file:
            rocket_inputs = json.load(file)
    except FileNotFoundError:
        return "Error: The specified JSON file does not exist."
    except json.JSONDecodeError:
        return "Error: The JSON file is not properly formatted."

    # validate required keys
    missing_keys = validate_required_parameters(rocket_inputs, required_keys)
    if missing_keys:
        return f"Missing required parameters: {', '.join(missing_keys)}"

    return rocket_inputs


# make sure all required parameters are present
def validate_required_parameters(parameters, required_keys):
    missing_keys = [key for key in required_keys if key not in parameters]
    return missing_keys  # Returns a list of missing keys, empty if all keys are present
