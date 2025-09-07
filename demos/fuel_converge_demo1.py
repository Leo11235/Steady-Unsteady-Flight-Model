'''
take some rocket inputs and convergence on a fuel mass using the steady-state algorithm
'''
import os, sys, json
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from steady_backend import steady_main
import output_handler

rocket_name = "Ancalagon_1"

simulation_type = {"simple fuel mass convergence": True}

simulation_settings = {
    "debug comments": "detailed comments", # "none", "comments", "detailed comments", or "very detailed comments"
    "output unit system": "MRT", # SI or MRT (rocket team mishmash of SI and IMP), or nice SI (some values get converted to cm or mm, etc;;; not yet created)
    "show graphs": True, # True or False, whether to show output graphs
    #"save simulation data": True, # currently only supports parametric tests
}

parametric_inputs={}

if __name__ == "__main__":
    input_file_path = f"./static/rocket_presets/{rocket_name}.json"
    with open(input_file_path, 'r') as json_file:
        rocket_inputs = json.load(json_file)
        
    output_dict = steady_main.run_steady(rocket_inputs, parametric_inputs, simulation_type, simulation_settings)
    output_handler.print_simulation_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict)