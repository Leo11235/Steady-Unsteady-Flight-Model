'''
create a 3-D graph of Isp vs ox flow rate on a range of 5-8 kg/s with a step size of 1, 
AND Fuel cell length from 40 to 100 cm with a step size of 10 cm
'''
import os, sys, json
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from steady_backend import steady_main
import output_handler

rocket_name = "Ancalagon_1"

simulation_type = {"parametric study": True}

simulation_settings = {
    "debug comments": "comments", # "none", "comments", "detailed comments", or "very detailed comments"
    "output unit system": "MRT", # SI or MRT (rocket team mishmash of SI and IMP), or nice SI (some values get converted to cm or mm, etc;;; not yet created)
    "show graphs": True, # True or False, whether to show output graphs
    "save simulation data": True, # currently only supports parametric tests
}

parametric_inputs = {
    "study name": "", # a name will be generated for any given parametric study, can also give a custom name
    "options of variables to parametrize": [ 
        # implemented, working
        "oxidizer mass flow rate", 
        "average oxidizer to fuel ratio", 
        "chamber pressure",
        "fuel length",
        "target apogee",
        "launch angle",
        "drag coefficient", 
        "dry mass",
        "rocket external diameter",
        "fuel grain density",
        
        # not yet implemented
        "nozzle throat radius", 
        "nozzle exit radius", 
        ],
    "variables to parametrize": [1,0,0,1,0,0,0,0,0], # either 1 or 0;;;; will eventually just change this to be directly a list of the actual variables we want to parametrize rather than just 1 or 0
    # variables below
    "oxidizer mass flow rate": { # 1
        "unit": "SI", # kg/s or lb/s
        "low end": 5,
        "high end": 8,
        "step size": 1,
    },
    "average oxidizer to fuel ratio": { # 2
        "unit": "dimensionless",
        "low end": 1,
        "high end": 5,
        "step size": 1,
    },
    "chamber pressure": { # 3
        "unit": "IMP", # Pa or psi
        "low end": 100,
        "high end": 2500,
        "step size": 100,
    },
    "fuel length": { # 4
        "unit": "SI", # m or in
        "low end": 40/100, # cm to m
        "high end": 100/100,
        "step size": 10/100,
    },
    "target apogee": { # 5
        "unit": "IMP", # m or ft
        "low end": 5000,
        "high end": 100000,
        "step size": 5000,
    },
    "launch angle": { # 6
        "unit": "dimensionless", # in degrees
        "low end": 2,
        "high end": 12,
        "step size": 2,
    },
    "drag coefficient": { # 7
        "unit": "dimensionless", # try values between .4 and 1
        "low end": .1,
        "high end": 3,
        "step size": .3,
    },
    "dry mass": { # 8
        "unit": "SI", # m or ft
        "low end": 50,
        "high end": 120,
        "step size": 10,
    },
    "rocket external diameter": { # 9
        "unit": "IMP", # m or in
        "low end": 6,
        "high end": 12,
        "step size": 2,
    },
    "fuel grain density": { # 10
        "unit": "SI", # m or ft
        "low end": 700,
        "high end": 2000,
        "step size": 200,
    },
    
    # below: not yet working
    "nozzle throat radius": {
        "unit": "SI", # m or in
        "low end": 10/1000, # mm to m
        "high end": 40/1000,
        "step size": 10/1000,
    },
    "nozzle exit radius": {
        "unit": "SI", # m or in
        "low end": 30/1000, # mm to m
        "high end": 50/1000,
        "step size": 10/1000,
    },
}

if __name__ == "__main__":
    input_file_path = f"./static/rocket_presets/{rocket_name}.json"
    with open(input_file_path, 'r') as json_file:
        rocket_inputs = json.load(json_file)
        
    output_dict = steady_main.run_steady(rocket_inputs, parametric_inputs, simulation_type, simulation_settings)
    output_handler.print_simulation_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict)