'''
this class handless the input file and format is in a way
that is nice to use for the rest of the program
'''

import json, re

simulation_settings_dict = {
    "simulation type" : "", 
    "debug comments": "normal comments", # "none", "normal comments", "detailed comments"
    "output unit system": "MRT", # SI or MRT (rocket team mishmash of SI and IMP)
    "show graphs": True, # whether to show output graphs (flight data, parametric study graphs, etc)
    "save simulation data": True, 
    "parametric study settings": {} # info gets filled in using input json
}

constants_dict = { # dictionary containing natural & simulation constants
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

    # simulation constants
    "number of timesteps": 100, # number of timesteps in burntime (only accounts for the part of the burn where there is a thrust). Ascent timesteps typically 300-400
    "tolerated apogee difference": .01, # m, allowed space between target and real apogee
    "smallest allowed inner fuel radius": 0.01, # m; equivalent to 1 cm
}

# process inputs, modifies simulation_settings_dict and constants_dict, and outputs a rocket_inputs dict
def read_input_file(input_file):
    # extract json from input file
    with open(input_file, 'r') as f:
        content = f.read()
    # remove comments
    content = re.sub(r'//.*?$|/\*.*?\*/', '', content, flags=re.MULTILINE | re.DOTALL)
    # parse the cleaned JSON string
    input_file =json.loads(content)
        
    # unpack input_file dicts
    simulation_settings_inputs = input_file['simulation settings']
    simulation_constants_inputs = input_file['simulation constants']
    rocket_inputs = input_file['rocket attributes']
        
    # override simulation settings and constants as needed        
    simulation_settings_dict.update({k: v for k, v in simulation_settings_inputs.items() if k in simulation_settings_dict})
    constants_dict.update({k: v for k, v in simulation_constants_inputs.items() if k in constants_dict})
    
    return rocket_inputs