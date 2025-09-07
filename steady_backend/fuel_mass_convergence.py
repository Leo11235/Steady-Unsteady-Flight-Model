# this file handles fuel mass convergence

from steady_backend.initialize_variables import error_dict
import steady_backend.calculations as calculations, steady_backend.pyPROPEP as runPROPEP, steady_backend.rocket_ascent_model as model_rocket_ascent

red = "\033[91m"     # Red text
reset_color = "\033[0m"    # Reset to default text color

# -------------------------------------------------------------
# -------------------- FUEL MASS ITERATION --------------------
# -------------------------------------------------------------

def iterate_over_inner_radius(rocket_input_values, constants_dict, simulation_settings):
    # create a dictionary which information about previously chosen inner radii and their apogees
    iteration_dict = {} # will add an entry after each iteration
    
    # set initial guess for internal radius (Ri0) using local variables
    lower_bound = 0 # lower bound for inner fuel radius
    upper_bound = rocket_input_values["fuel external diameter"] / 2 # diameter -> radius; upper bound for inner fuel radius
    Ri0 = (upper_bound + lower_bound) / 2 # set an initial guess, from testing it's not worth bothering to make a better guess
    
    if constants_dict["sea level gravity"] < 1: # if gravity is very low, get a very low fuel mass
        # this is almost never true, but necessary for some of the meme planets (ie death star and B612)
        Ri0 = rocket_input_values["fuel external diameter"] / 2 - constants_dict["sea level gravity"] / 100
        
    # iteration goes here
    i = 1 # num of iterations
    rocket_parameters = { 'reached apogee' : 0 } # bogus value so the while loop can start
    while(abs(rocket_parameters["reached apogee"] - rocket_input_values["target apogee"]) >= constants_dict["tolerated apogee difference"]):
        # clear rocket_parameters, add Ri0 guess
        rocket_parameters = { "initial internal fuel radius": Ri0 }
        # run calculations
        calculations.CV2_calculations(rocket_input_values, rocket_parameters)
        runPROPEP.run_PROPEP(rocket_input_values, rocket_parameters)
        calculations.CV3_calculations(rocket_input_values, rocket_parameters, constants_dict)
        # model rocket ascent
        flight_dict = model_rocket_ascent.simulate_rocket_ascent(rocket_input_values, rocket_parameters, constants_dict) # create a dictionary with tons of information on forces, kinematics, etc with respect to time
        # refine Ri0 guess
        rocket_parameters["reached apogee"] = flight_dict["altitude"][-1]
        if rocket_parameters["reached apogee"] > rocket_input_values["target apogee"]:
            lower_bound = Ri0 # if the rocket flies too high, increase inner radius (less fuel)
        elif rocket_parameters["reached apogee"] < rocket_input_values["target apogee"]:
            upper_bound = Ri0 # & vice versa
        Ri0 = (upper_bound + lower_bound) / 2
        
        # add to iteration dict
        rocket_parameters["reached apogee"] = 0 if rocket_parameters["reached apogee"] < 0 else rocket_parameters["reached apogee"] # make the rocket not phase through the ground
        apogee_diff = abs(rocket_parameters["reached apogee"] - rocket_input_values["target apogee"])
        relative_apogee_diff = apogee_diff / rocket_input_values["target apogee"]
        iteration_dict[f'Iteration {i}'] = f'{value_to_rgb(relative_apogee_diff)}Fuel inner radius: {round(rocket_parameters["initial internal fuel radius"] * 39.3701, 3):g} (in) --> apogee: {round(rocket_parameters["reached apogee"] * 3.28084, 1):g} (ft){reset_color}'
        print(f"  {i}  " + iteration_dict[f'Iteration {i}'] + (f" ... {red}rocket failed to launch!{reset_color}" if rocket_parameters["reached apogee"] == 0 else "")) if simulation_settings["debug comments"] == "detailed comments" else None
                
        i+=1
        
        # if the fuel mass is as high as it can get and the rocket's reached apogee is significantly lower than target apogee, return error
        if rocket_parameters["initial internal fuel radius"] == constants_dict["smallest allowed inner fuel radius"] and rocket_input_values["target apogee"] > rocket_parameters["reached apogee"] + 1000:
            error_dict["error"] = True
            error_dict["bisection error"] = True
            error_dict["bisection error information"] = f"During the fuel convergence process, the rocket's inner fuel cell radius (currently {round(rocket_parameters["initial internal fuel radius"], 2)*100} cm) has fallen to a value below {constants_dict["smallest allowed inner fuel radius"]*100} cm. This indicates that either the rocket's input values or the target apogee are unrealistic expectations for this given set of inputs."
            print("\033[38;2;255;165;0mFuel convergence error type 1\033[0m") if simulation_settings["debug comments"] == "detailed comments" else None
            print(error_dict["bisection error information"]) if simulation_settings["debug comments"] == "very detailed comments" else None
            return fuel_convergence_success(False, rocket_input_values, constants_dict, rocket_parameters, flight_dict)
        
        # quick check to make sure the inner fuel radius is not too small (& therefore the rocket has no hope of ever reaching apogee)
        if rocket_parameters["initial internal fuel radius"] < constants_dict["smallest allowed inner fuel radius"]:
            Ri0 = constants_dict["smallest allowed inner fuel radius"] # MIGHT CAUSE PROBLEM
            print("\033[38;2;255;165;0mFuel convergence error type 2 (inner radius dipped below min allowed value, no cause for concern if error type 1 appears after next convergence)\033[0m") if simulation_settings["debug comments"] == "detailed comments" else None
            #error_dict["error"] = True
            #error_dict["bisection error"] = True
            #error_dict["bisection error information"] = f"During the fuel convergence process, the rocket's inner fuel cell radius (currently {round(rocket_parameters["initial internal fuel radius"], 2)*100} cm) has fallen to a value below {constants_dict["smallest allowed inner fuel radius"]*100} cm. This indicates that either the rocket's input values or the target apogee are unrealistic expectations for this given set of inputs."
            #print(error_dict["bisection error information"]) if simulation_settings["debug comments"] == "comments" else None
            #rocket_parameters["apogee reached T/F"] = False
            #return {"rocket inputs": rocket_input_values, 
            #        "constants dict": constants_dict, 
            #        "rocket parameters": rocket_parameters, 
            #        "flight data": flight_dict}
            
        # exit the program in cases where the inner fuel radius keeps failing to converge
        if i > 30 and rocket_parameters["reached apogee"] < rocket_input_values["target apogee"]:
            print("\033[91mFuel convergence error type 3 (this one is kinda bad and should hopefully never happen but also is not too terrible if you are running parametric studies and most of the 'reached target apogee' statuses are false before this but true after (in this case the rocket is almost able to reach the target apogee but not quite there))\033[0m")
            return fuel_convergence_success(False, rocket_input_values, constants_dict, rocket_parameters, flight_dict)
        
    # apogee has been reached at this point
    return fuel_convergence_success(True, rocket_input_values, constants_dict, rocket_parameters, flight_dict)
    
# helper function to return end value from fuel convergence algorithm, exists because there are several exit conditions for the function above
def fuel_convergence_success(convergence_successful: bool, rocket_input_values, constants_dict, rocket_parameters, flight_dict):
    if convergence_successful:
        rocket_parameters["apogee reached T/F"] = True
    else: # for whatever reason, the fuel convergence failed (usually bissection error)
        rocket_parameters["apogee reached T/F"] = False
    return {"rocket inputs": rocket_input_values, 
            "constants dict": constants_dict, 
            "rocket parameters": rocket_parameters, 
            "flight data": flight_dict}


# helper functions

# takes a float (abs value difference between target apogee and real achieved apogee, divided by target apogee), returns an RGB color
def value_to_rgb(value):
    if value > 0.1:
        rgb = (255, 0, 0)  # Very Red
    elif value > 0.01:
        rgb = (255, 69, 0)  # Red-Orange
    elif value > 0.001:
        rgb = (255, 165, 0)  # Orange
    elif value > 0.0001:
        rgb = (255, 255, 0)  # Yellow
    elif value > 0.00001:
        rgb = (173, 255, 47)  # Yellow-Green
    elif value > 0.000001:
        rgb = (50, 205, 50)  # Green
    else:
        rgb = (0, 255, 0)  # Very Green
    return f"\033[38;2;{rgb[0]};{rgb[1]};{rgb[2]}m"