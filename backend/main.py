'''
guide for dictionaries: 
    rocket_inputs refers to the user's inputs, these will never change between each simulation
    rocket_parameters refers to the parameters that get calculated in each iteration of the simulation, only the final iteration's parameters are kept
    constants_dict contains universal values such as sea level gravity
    flight_dict contains n-long lists, where n is the amount of timesteps; these will house information for graph creation
    burn_dict, ascent_dict respectively contain information about the rocket's flight during the burn and after the burn (in file rocket_ascent_simulator.py)
'''

from variable_initialization import simulation_settings_dict, read_input_file
from simulation_iterations import iterate_over_fuel_mass, simulate_rocket_burn

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
        print(generate_parametric_study_description())
        return
    
    # optimize values for unsteady
    elif simulation_settings_dict["simulation type"].lower() == "optimize values for unsteady": # optimizes Isp as a function of Mo, Lf, Re for use by unsteady -- basically a specific parametric study
        return


def generate_parametric_study_description():
    # construct a phrase for each parameter
    phrases = []
    for key, value in simulation_settings_dict["parametric study settings"].items():
        phrase = f"\n   {key} (low end: {value['low end']}, high end: {value['high end']}, step size: {value['step size']})"
        phrases.append(phrase)
    
    # join with commas and 'and' before the last item
    if len(phrases) == 1:
        summary = phrases[0]
    elif len(phrases) == 2:
        summary = " and ".join(phrases)
    else:
        summary = ", ".join(phrases[:-1]) + ", and " + phrases[-1]
    
    return f"Running parametric study on {summary}"
    

# prints all items of a dictionary in a nice way
def print_dict(input_dictionary):
    print()
    max_key_length = max(len(str(key)) for key in input_dictionary.keys())
    for key, value in input_dictionary.items():
        print(f"{key} {"." * (max_key_length - len(str(key)) + 2)} {value}")
    print()

if __name__ == '__main__':
    file = "./input_files/Esteban's_Ancalagon_param_1.jsonc"
    main(file)
