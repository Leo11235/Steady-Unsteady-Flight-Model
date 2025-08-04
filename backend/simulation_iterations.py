from variable_initialization import simulation_settings_dict, constants_dict
import calculations
from pyPROPEP import runPROPEP
from rocket_ascent_simulator import simulate_rocket_ascent

# converges on an ideal fuel mass given a target apogee; returns rocket_inputs, rocket_parameters, flight_dict
def iterate_over_fuel_mass(rocket_inputs):
    # set initial bogus values so the while loop can start
    lower_bound = 0 # lower bound for inner fuel radius
    upper_bound = rocket_inputs["fuel external radius"] # upper bound
    rocket_parameters = {
        "reached apogee" : float("-inf"), # bogus value so the while loop can start
        "initial internal fuel radius": rocket_inputs["fuel external radius"]/2 # initial internal fuel radius guess     
    }
    
    # run simulation with a certain fuel mass, check apogee, refine guess, run again, .... until the final apogee reaches the desired value
    i=1
    while(True):
        rocket_parameters = simulate_rocket_burn(rocket_inputs, rocket_parameters)
        flight_dict = simulate_rocket_ascent(rocket_inputs, rocket_parameters, constants_dict)
        
        rocket_parameters["reached apogee"] = flight_dict["altitude"][-1]
        if correct_apogee_reached(rocket_inputs, rocket_parameters):
            print(f'SIMULATION TERMINATED AFTER {i} LOOPS. Apogee achieved: {round(rocket_parameters["reached apogee"], 1)} meters ({round(abs(rocket_parameters["reached apogee"] - rocket_inputs["target apogee"]), 5)} meters off from target)') if "comments" in simulation_settings_dict["debug comments"] else None
            return rocket_inputs, rocket_parameters, flight_dict
        else: # refine inner fuel radius guess and try again
            if rocket_parameters["reached apogee"] > rocket_inputs["target apogee"]:
                lower_bound = rocket_parameters["initial internal fuel radius"] # if the rocket flies too high, increase inner radius (less fuel)
            elif rocket_parameters["reached apogee"] < rocket_inputs["target apogee"]:
                upper_bound = rocket_parameters["initial internal fuel radius"] # & vice versa
            else: 
                print("ERROR: this should be impossible")
        
        print(f'Loop {i} of simulation, apogee achieved: {round(rocket_parameters["reached apogee"])} meters') if simulation_settings_dict["debug comments"] == "detailed comments" else None
        
        # reset rocket parameters dict before re-launching simulation
        rocket_parameters = {"initial internal fuel radius": (upper_bound + lower_bound) / 2}
        i+=1

# helper function to check whether the desired apogee was reached for the rocket
def correct_apogee_reached(rocket_inputs, rocket_parameters):
    return True if abs(rocket_parameters["reached apogee"] - rocket_inputs["target apogee"]) <= constants_dict["tolerated apogee difference"] else None



# takes rocket input data, outputs rocket performance becnhmarks
def simulate_rocket_burn(rocket_inputs, rocket_parameters):
    # step 1
    calculations.CV2_calculations(rocket_inputs, rocket_parameters)
    # step 2: propep
    runPROPEP(rocket_inputs, rocket_parameters)
    # step 3
    calculations.CV3_calculations(rocket_inputs, rocket_parameters, constants_dict)
    # output
    return rocket_parameters

