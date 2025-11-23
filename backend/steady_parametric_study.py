import itertools, numpy as np
from . import steady_simulation_iterations as steady_simulation_iterations
from .steady_variable_initialization import simulation_settings_dict

# runs parametric study, works for n inputs
def run_paramatric_study(rocket_inputs, simulation_settings_dict):
    param_settings = simulation_settings_dict["parametric study settings"]
    
    # generate ranges for each variable
    var_ranges = {}
    for var_name, var_values in param_settings.items():
        var_ranges[var_name] = generate_range(var_values["low end"], var_values["high end"], var_values["step size"], simulation_settings_dict.get("parametric study step size undershoot tolerance", 0.01))
        
    # calculate total number of iterations
    total_iterations = np.prod([len(r) for r in var_ranges.values()])
    print(f"Total combinations to run: {total_iterations}\n")    
    
    # get variable names
    variables = list(param_settings.keys())
    # create dict to hold results
    param_results_dict = {
        "variable ranges": var_ranges, 
        "rocket inputs": [], 
        "rocket parameters": [], 
        "flight dict": []
        }
    
    # iterate through all combinations
    i = 0
    for combination in itertools.product(*var_ranges.values()):
        i+=1
        # update rocket inputs
        current_rocket_inputs = rocket_inputs.copy()
        for var, val in zip(variables, combination):
            current_rocket_inputs[var] = val
        # run simulation
        print(f'{round(100*i/total_iterations, 1)}% -- Loop {i} {combination}') if "comments" in simulation_settings_dict["debug comments"] else None
        rocket_parameters, flight_dict = steady_simulation_iterations.iterate_over_fuel_mass(current_rocket_inputs)
        # save results
        param_results_dict["rocket inputs"].append(current_rocket_inputs)
        param_results_dict["rocket parameters"].append(rocket_parameters)
        param_results_dict["flight dict"].append(flight_dict)
        
        # STILL LEFT TO ADD
        # save every 100 iterations
    
    return param_results_dict

# helper function that generates a statement describing the parametric study
def generate_parametric_study_description(simulation_settings_dict):
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
    
    return f"Running parametric study on {summary}\n"

# helper function to generate range of values to iterate over
def generate_range(low, high, step, tolerance):
    values = [low]
    while True:
        next_value = values[-1] + step
        if next_value > high: # if next value exceeds high
            # check if current value is close enough to high
            if (high - values[-1]) < tolerance * step:
                break  # close enough, don't add next_value
            else:
                values.append(next_value)  # overshoot is acceptable
                break
        else:
            values.append(next_value)
    return values

