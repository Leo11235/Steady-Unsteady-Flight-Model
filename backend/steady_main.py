'''
guide for dictionaries: 
    rocket_inputs refers to the user's inputs, these will never change between each simulation
    rocket_parameters refers to the parameters that get calculated in each iteration of the simulation, only the final iteration's parameters are kept
    constants_dict contains universal values such as sea level gravity
    flight_dict contains n-long lists, where n is the amount of timesteps; these will house information for graph creation
    burn_dict, ascent_dict respectively contain information about the rocket's flight during the burn and after the burn (in file rocket_ascent_simulator.py)
'''

from steady_variable_initialization import simulation_settings_dict, read_input_file
from steady_simulation_iterations import iterate_over_fuel_mass, simulate_rocket_burn
import steady_parametric_study as steady_parametric_study
import json

def main(input_file_path):
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
        return simulate_rocket_burn(rocket_inputs, rocket_parameters)
    
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

if __name__ == '__main__':
    # get file and run simulation
    file = "./steady_input_files/Esteban's_Ancalagon.jsonc"
    data = main(file)
    
    # save data
    if simulation_settings_dict["save simulation data"]:
        with open("steady_output_files/simulation_output.json", "w") as f:
            json.dump(data, f, indent=4)
            
    # graph data
    
    # estimate 'file size' of data
    json_str = json.dumps(data)
    size_bytes = len(json_str.encode('utf-8'))
    size_kb = size_bytes / 1024
    size_mb = size_kb / 1024
    print(f"Estimated JSON size: {size_bytes} bytes ({size_kb:.2f} KB / {size_mb:.2f} MB)")

