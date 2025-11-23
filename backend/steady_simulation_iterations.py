from .steady_variable_initialization import simulation_settings_dict, constants_dict
from . import steady_calculations as steady_calculations
from .pyPROPEP import runPROPEP
from .steady_rocket_ascent_simulator import simulate_rocket_ascent

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
    smallest_allowed_fuel_radius_tested = False # useful for when the rocket is just barely able to scrape it to apogee
    while(True):
        # simulate rocket
        rocket_parameters = simulate_rocket_burn(rocket_inputs, rocket_parameters)
        flight_dict = simulate_rocket_ascent(rocket_inputs, rocket_parameters, constants_dict)
        # get apogee
        rocket_parameters["reached apogee"] = flight_dict["altitude"][-1]
        if correct_apogee_reached(rocket_inputs, rocket_parameters):
            print(f'    LOOP {i} - FINAL APOGEE: {round(rocket_parameters["reached apogee"], 1)} meters ({round(abs(rocket_parameters["reached apogee"] - rocket_inputs["target apogee"]), 5)} meters off from target)') if "detailed comments" in simulation_settings_dict["debug comments"] else None
            return {"rocket parameters": rocket_parameters, "flight dict": flight_dict}
        else: 
            # check if rocket failed to reach apogee even with lowest allowed internal fuel radius
            if rocket_parameters["initial internal fuel radius"] == constants_dict["smallest allowed inner fuel radius"]: # AND not reached apogee, but that possibility has already been filtered out
                # if the code makes it to here, the rocket cannot reach the target apogee
                print(f"    ROCKET CANNOT REACH APOGEE (final apogee reached: {rocket_parameters["reached apogee"]} meters)")
                return {"rocket parameters": rocket_parameters, "flight dict": flight_dict}
            # refine inner fuel radius guess and try again
            if rocket_parameters["reached apogee"] > rocket_inputs["target apogee"]:
                lower_bound = rocket_parameters["initial internal fuel radius"] # if the rocket flies too high, increase inner radius (less fuel)
            elif rocket_parameters["reached apogee"] < rocket_inputs["target apogee"]:
                upper_bound = rocket_parameters["initial internal fuel radius"] # & vice versa
            else: 
                print("ERROR: this should be impossible")
        
        print(f'    Loop {i} of simulation, apogee achieved: {round(rocket_parameters["reached apogee"], 1)} meters') if simulation_settings_dict["debug comments"] == "detailed comments" else None
        
        # reset rocket parameters and re-guess initial internal fuel radius
        rocket_parameters = {"initial internal fuel radius": (upper_bound + lower_bound) / 2}
        
        # if the initial internal fuel radius dips below the smallest allowed amount, run the simulation with the smallest allowed amount instead
            # if that still fails, that means the rocket cannot reach the apogee and the simulation is aborted
        if rocket_parameters["initial internal fuel radius"] < constants_dict["smallest allowed inner fuel radius"]:
            rocket_parameters["initial internal fuel radius"] = constants_dict["smallest allowed inner fuel radius"]
            lower_bound = constants_dict["smallest allowed inner fuel radius"]
        
        i+=1

# helper function to check whether the desired apogee was reached for the rocket
def correct_apogee_reached(rocket_inputs, rocket_parameters):
    return True if abs(rocket_parameters["reached apogee"] - rocket_inputs["target apogee"]) <= constants_dict["tolerated apogee difference"] else None

# takes rocket input data, outputs rocket performance becnhmarks
# used mainly for hotfire tests
def simulate_rocket_burn(rocket_inputs, rocket_parameters):
    # step 1
    steady_calculations.CV2_calculations(rocket_inputs, rocket_parameters)
    # step 2: propep
    runPROPEP(rocket_inputs, rocket_parameters)
    # step 3
    steady_calculations.CV3_calculations(rocket_inputs, rocket_parameters, constants_dict)
    # output
    return rocket_parameters