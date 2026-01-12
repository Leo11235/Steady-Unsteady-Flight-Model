'''
guide for dictionaries: 
    rocket_inputs refers to the user's inputs, these will never change between each simulation
    rocket_parameters refers to the parameters that get calculated in each iteration of the simulation, only the final iteration's parameters are kept
    constants_dict contains universal values such as sea level gravity
    flight_dict contains n-long lists, where n is the amount of timesteps; these will house information for graph creation
    burn_dict, ascent_dict respectively contain information about the rocket's flight during the burn and after the burn (in file rocket_ascent_simulator.py)
'''

from .steady_variable_initialization import simulation_settings_dict, read_input_file
from .steady_simulation_iterations import iterate_over_fuel_mass, simulate_rocket_burn
from . import steady_parametric_study as steady_parametric_study
import json

def steady_main(input_file_path):
    # create dictionaries of rocket inputs, simulation settings
    rocket_inputs = read_input_file(input_file_path)
        
    # determine simulation type and run the appropriate one 
    # hotfire test
    if simulation_settings_dict["simulation type"].lower() == "hotfire": # calculates internal rocket dynamics without calculating rocket flight
        print('Running hotfire test')
        if not "initial internal fuel radius" in rocket_inputs or not rocket_inputs["initial internal fuel radius"]: # requires 'initial internal fuel radius' to be in rocket_inputs
            print("ERROR: Must enter a valid initial internal fuel radius!")
            exit()
        # put initial internal fuel radius into rocket_parameters for later calculations
        rocket_parameters = {"initial internal fuel radius": rocket_inputs["initial internal fuel radius"]} 
        # run simulation
        return rocket_inputs, simulate_rocket_burn(rocket_inputs, rocket_parameters)
    
    # fuel mass convergence
    elif simulation_settings_dict["simulation type"].lower() == "fuel mass convergence": # converges on ideal fuel mass given a target apogee
        print("Running fuel mass convergence simulation")
        return iterate_over_fuel_mass(rocket_inputs)
    
    # parametric study
    elif simulation_settings_dict["simulation type"].lower() == "parametric study": # runs many fuel mass convergences
        print(steady_parametric_study.generate_parametric_study_description(simulation_settings_dict))
        return steady_parametric_study.run_paramatric_study(rocket_inputs, simulation_settings_dict)
        
    
    # optimize values for unsteady
    elif simulation_settings_dict["simulation type"].lower() == "optimize values for unsteady": # optimizes Isp as a function of Mo, Lf, Re for use by unsteady -- basically a specific parametric study
        return
    

# prints all items of a dictionary in a nice way
def print_dict(input_dictionary):
    print()
    max_key_length = max(len(str(key)) for key in input_dictionary.keys())
    for key, value in input_dictionary.items():
        print(f"{key} {"." * (max_key_length - len(str(key)) + 2)} {value}")
    print()
    
def format_hotfire_results(rocket_inputs, rocket_parameters):
    hotfire_dict = {
        "rocket name": rocket_inputs["rocket name"], 
        "thrust": rocket_parameters["thrust"], 
        "burntime": rocket_parameters["burntime"],
        "Isp": rocket_parameters["Isp"], 
        "total impulse": rocket_parameters["total impulse"], 
        "average oxidizer mass flow rate": rocket_inputs["oxidizer mass flow rate"],
        "chamber pressure": rocket_inputs["chamber pressure"],
        "oxidizer mass": rocket_parameters["burntime"] * rocket_inputs["oxidizer mass flow rate"],
        "fuel mass": rocket_parameters["fuel mass"],
        "fuel length": rocket_inputs["fuel length"],
        "fuel external radius": rocket_inputs["fuel external radius"], 
        "fuel internal radius": rocket_parameters["initial internal fuel radius"],
        "chamber temperature": rocket_parameters["chamber temperature"], 
        
    }
    print("=================================================\n========== STEADY HOTFIRE TEST RESULTS ==========\n=================================================\n")
    print_dict(hotfire_dict)
    
    return hotfire_dict
