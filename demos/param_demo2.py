'''
create a 2-D graph of Isp vs chamber pressure on a range of 300-700 psi with a step size of 100
'''
import os, sys, json
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from steady_backend import steady_main
import output_handler

rocket_name = "Ancalagon_1"

simulation_type = {"parametric study": True}

simulation_settings = {
    "debug comments": "detailed comments", # "none", "comments", "detailed comments", or "very detailed comments"
    "output unit system": "MRT", # SI or MRT (rocket team mishmash of SI and IMP), or nice SI (some values get converted to cm or mm, etc;;; not yet created)
    "show graphs": True, # True or False, whether to show output graphs
    "save simulation data": False, # currently only supports parametric tests
}

parametric_inputs = {
    "options of variables to parametrize": [ 
        # implemented, working
        "oxidizer mass flow rate", 
        "average oxidizer to fuel ratio", 
        "chamber pressure",
        ],
    "variables to parametrize": [0,0,1,0,0,0,0,0,0], # either 1 or 0;;;; will eventually just change this to be directly a list of the actual variables we want to parametrize rather than just 1 or 0
    # variables below
    "oxidizer mass flow rate": { # 1
        "unit": "SI", # kg/s or lb/s
        "low end": 1,
        "high end": 80,
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
        "low end": 300,
        "high end": 700,
        "step size": 100,}}

if __name__ == "__main__":
    input_file_path = f"./static/rocket_presets/{rocket_name}.json"
    with open(input_file_path, 'r') as json_file:
        rocket_inputs = json.load(json_file)
        
    output_dict = steady_main.run_steady(rocket_inputs, parametric_inputs, simulation_type, simulation_settings)
    output_handler.print_simulation_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict)