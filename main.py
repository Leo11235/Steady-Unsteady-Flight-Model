'''
guide for dictionaries: 
    rocket_inputs refers to the user's inputs, these will never change between each simulation
    rocket_parameters refers to the parameters that get calculated in each iteration of the simulation, only the final iteration's parameters are kept
    constants_dict contains universal values such as sea level gravity
    timestep_dict contains n-long lists, where n is the amount of timesteps; these will house information for graph creation
'''

from variable_initialization import constants_dict, read_rocket_attributes
import calculations
from pyPROPEP import runPROPEP
import rocket_ascent_model

def main():
    # create dictionary of rocket inputs
    rocket_inputs = read_rocket_attributes("./input_files/Esteban's_ancalagon.json")    
    # create a dictionary of rocket parameters, these values will change after every iteration
    rocket_parameters = {}

    # iterate over inner radius

    # TESTING: one single iteration
    rocket_parameters["initial internal fuel radius"] = 0.039974 # random inner radius for testing
    # step 1
    calculations.CV2_calculations(rocket_inputs, rocket_parameters) # run CV2 with this r
    # step 2
    runPROPEP(rocket_inputs, rocket_parameters)
    # step 3
    calculations.CV3_calculations(rocket_inputs, rocket_parameters, constants_dict)
    # step 4: CV4 calculations
    rocket = rocket_ascent_model.Rocket(rocket_inputs, rocket_parameters)


    print_dict(rocket_inputs)
    print_dict(rocket_parameters)



    # altitude convergence (iteration)
        # set bounds for initial internal fuel radius
        # calculations part 1
        # PROPEP
        # calculations part 2
        # time functions (more calculations)
    
    # parametric studies

    # put relevant results in a JSON file



# prints all items of a dictionary in a nice way
def print_dict(input_dictionary):
    print()
    max_key_length = max(len(str(key)) for key in input_dictionary.keys())
    for key, value in input_dictionary.items():
        fancyshmancy = "." * (max_key_length - len(str(key)) + 2) 
        print(f"{key} {fancyshmancy} {value}")
    print()

if __name__ == '__main__':
    main()