'''
guide for dictionaries: 
    rocket_inputs refers to the user's inputs, these will never change between each simulation
    rocket_parameters refers to the parameters that get calculated in each iteration of the simulation, only the final iteration's parameters are kept
    constants_dict contains universal values such as sea level gravity
    overall_flight_dict contains n-long lists, where n is the amount of timesteps; these will house information for graph creation
    burn_dict, ascent_dict respectively contain information about the rocket's flight during the burn and after the burn
'''

from variable_initialization import constants_dict, read_rocket_attributes
import calculations
from pyPROPEP import runPROPEP
import rocket_ascent_model
from graphing import graph_everything

def main():
    # create dictionary of rocket inputs
    rocket_inputs = read_rocket_attributes("./input_files/Esteban's_ancalagon.json")    

    # iterate over fuel inner radius to eventually reach the right apogee
    rocket_parameters, overall_flight_dict = iterate_over_fuel_mass(rocket_inputs, constants_dict)

    # at this point, rocket_inputs and rocket_parameters contain all the information we need about the rocket, 
    # and overall_flight_dict contains info that can be used to graph forces and kinematics

# still need to implement:
    # parametric studies
    # put relevant results in a JSON file


    
    # DEBUGGING 
    graph_everything(overall_flight_dict)
    print_dict(rocket_inputs)
    print_dict(rocket_parameters)
    #print_dict(timesteps_dict)
    #print_dict(ascent_dict)



# ITERATION
# inner radius determines fuel mass -> fuel mass determines rocket's performance -> the rocket reaches a certain altitude
def iterate_over_fuel_mass(rocket_inputs, constants_dict):
    # create a dictionary of rocket parameters, these values will change after every iteration
    rocket_parameters = {}

    # set initial guess for internal radius using local variables
    lower_bound = 0
    upper_bound = rocket_inputs["fuel external diameter"] / 2 # diameter -> radius
    Ri0 = (upper_bound + lower_bound) / 2

    # iteration goes here
    i = 1 # num of iterations
    rocket_parameters = { 'reached apogee' : 0 } # bogus value so the while loop can start
    while(abs(rocket_parameters["reached apogee"] - rocket_inputs["target apogee"]) >= constants_dict["tolerated apogee difference"]):
        # step 0: clear dictionaries, set lower bound
        rocket_parameters = { "initial internal fuel radius": Ri0 }        
        # step 1
        calculations.CV2_calculations(rocket_inputs, rocket_parameters) # run CV2 with this r
        # step 2: propep
        runPROPEP(rocket_inputs, rocket_parameters)
        # step 3
        calculations.CV3_calculations(rocket_inputs, rocket_parameters, constants_dict)
        # step 4: CV4 calculations
        rocket = rocket_ascent_model.Rocket(rocket_inputs, rocket_parameters) # create a rocket object with all the atributes previously calculated
        burn_dict, ascent_dict = rocket.calculate_kinematics(rocket_inputs, rocket_parameters, constants_dict) # calculate burn and ascent
        overall_flight_dict = rocket_ascent_model.merge_flight_dictionaries(burn_dict, ascent_dict) # create a dictionary with tons of information on forces, kinematics, with respect to time
        # step 5: refine inner radius guess; adjust bounds
        if rocket_parameters["reached apogee"] > rocket_inputs["target apogee"]:
            lower_bound = Ri0 # if the rocket flies too high, increase inner radius (less fuel)
        elif rocket_parameters["reached apogee"] < rocket_inputs["target apogee"]:
            upper_bound = Ri0 # & vice versa
        Ri0 = (upper_bound + lower_bound) / 2
        # quick check to make sure the inner fuel radius is not too small (& therefore the rocket has no hope of ever reaching apogee)
        if rocket_parameters["initial internal fuel radius"] < constants_dict["smallest allowed inner fuel radius"]:
            erocketyle_dysfunction() # exits the program
        i += 1
    return rocket_parameters, overall_flight_dict


# for when the rocket just can't make it :(
def erocketyle_dysfunction():
    print(
        "Rocket simulation halted: The rocket's parameters are unrealistic. "
        "The initial internal fuel radius has converged to a value smaller than the allowed minimum, "
        "indicating the rocket cannot reach the target apogee with the current setup."
    )
    # finish
    exit(0)


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
