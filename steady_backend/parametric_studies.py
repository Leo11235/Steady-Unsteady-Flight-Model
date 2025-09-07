import steady_backend.fuel_mass_convergence as fuel_mass_convergence, steady_backend.unit_conversion as convert
import numpy, itertools, inspect

# ANSI color codes for printing status updates
green = "\033[92m"   # Green text
red = "\033[91m"     # Red text
yellow = "\033[93m"  # Yellow text
blue = "\033[94m"    # Blue text
orange = "\033[38;5;214m"  # Orange text (approximation using 256-color mode)
reset = "\033[0m"    # Reset to default text color

def run_parametric_study(rocket_input_values, constants_dict, parametric_inputs, simulation_settings):
    # convert input values to SI
    convert.convert_param_to_SI(parametric_inputs)
    
    # run study
    variables_being_parametrized = sum(parametric_inputs["variables to parametrize"])
    if variables_being_parametrized == 0:
        print('Please input a variable to parametrize over')
        exit()
    elif variables_being_parametrized == 1:
        parametric_results_dict = run_param_one_input(rocket_input_values, constants_dict, parametric_inputs, simulation_settings)
    elif variables_being_parametrized == 2:
        parametric_results_dict = run_param_two_inputs(rocket_input_values, constants_dict, parametric_inputs, simulation_settings)
    elif variables_being_parametrized == 3:
        parametric_results_dict = run_param_three_inputs(rocket_input_values, constants_dict, parametric_inputs, simulation_settings)
    else:
        print(f"\033[91mParametric study error type 1 in {inspect.getfile(inspect.currentframe())}, Line: {inspect.currentframe().f_lineno} (bad)\033[0m")
        return None

    return parametric_results_dict

# -----------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY 1: One Input --------------------
# -----------------------------------------------------------------------

def run_param_one_input(rocket_input_values, constants_dict, parametric_inputs, simulation_settings):
    # determine variable to parametrize
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    variable_to_parametrize = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]]
    # create dictionary with bounds & step size
    param_dict = {
        "study name": parametric_inputs.get("study name"),
        f'{variable_to_parametrize} unit': parametric_inputs[variable_to_parametrize]["unit"],
        f'{variable_to_parametrize} low end': parametric_inputs[variable_to_parametrize]["low end"],
        f'{variable_to_parametrize} high end': parametric_inputs[variable_to_parametrize]["high end"],
        f'{variable_to_parametrize} step size': parametric_inputs[variable_to_parametrize]["step size"],
    }
    # create dictionary which will hold all the data points for this given parametric study
    param_results_dict = {
        f'{variable_to_parametrize}': [], # x-axis
        'specific impulse': [], # y-axis
        'apogee reached': [], # a list of T/F values indicating whether the target apogee was reached with this input
        'rocket parameters for given input': [] # a list of dicts, containing details about the rocket for each iteration
    }
    # Mo and OF are mututally exclusive: if OF is selected, must remove Mo from rocket inputs
    remove_unnecessary_rocket_input_vars(rocket_input_values, variable_to_parametrize)
    # generate list of input values to traverse over
    x_vals = numpy.arange(param_dict[f'{variable_to_parametrize} low end'], param_dict[f'{variable_to_parametrize} high end'] + param_dict[f'{variable_to_parametrize} step size'], param_dict[f'{variable_to_parametrize} step size'])
    # estimate runtime
    print(f"{blue}Performing parametric analysis on {variable_to_parametrize}: need to converge on target apogee {len(x_vals)} times{reset}") if simulation_settings["debug comments"] else None
    i=1
    for x in x_vals:
        # initialize rocket inputs dict
        rocket_input_values[variable_to_parametrize] = x
        # run fuel mass convergence algo
        fuel_convergence_results = fuel_mass_convergence.iterate_over_inner_radius(rocket_input_values, constants_dict, simulation_settings)
        # get Isp & rocket parameters
        rocket_parameters = fuel_convergence_results["rocket parameters"]
        Isp = rocket_parameters["Isp"]
        # add new values to param_results_dict
        param_results_dict[f'{variable_to_parametrize}'].append(x)
        param_results_dict["specific impulse"].append(Isp)
        param_results_dict["apogee reached"].append(rocket_parameters["apogee reached T/F"])
        param_results_dict["rocket parameters for given input"].append(rocket_parameters)
        
        # print status 
        apogee_reached = rocket_parameters["apogee reached T/F"]
        color = green if apogee_reached else red
        print(f'{round(100 * i/(len(x_vals)),1)}% -- {blue}{variable_to_parametrize}: {x}, Isp: {Isp}, reached apogee target: {color}{apogee_reached}{reset}') if simulation_settings["debug comments"] != "none" else None
        
        # save values and print status report if i is divisible by 100 (every 100 iterations)
        #if i % 100 == 0:
        #    print_param_hundredly_status_report(i, param_results_dict, round(100 * i/(len(x_vals)),1))
        #    output_handler.save_in_progress_param_simulation(param_dict, param_results_dict)
        
        i+=1
        
    return {
        "parametric study input dict": param_dict, 
        "parametric study output dict": param_results_dict}
    
    
# ------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY 2: Two Inputs --------------------
# ------------------------------------------------------------------------

def run_param_two_inputs(rocket_input_values, constants_dict, parametric_inputs, simulation_settings):
    # determine variable to parametrize
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    first_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]]
    second_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[1]]
    # check that input variables are not mutually exclusive
    param_study_inputs_not_mutually_exclusive_checker(first_var, second_var)
    # create dictionary with bounds & step size
    param_iteration_dict = {
        "study name": parametric_inputs.get("study name"),
        f'{first_var} unit': parametric_inputs[first_var]["unit"],
        f'{first_var} low end': parametric_inputs[first_var]["low end"],
        f'{first_var} high end': parametric_inputs[first_var]["high end"],
        f'{first_var} step size': parametric_inputs[first_var]["step size"],
        f'{second_var} unit': parametric_inputs[second_var]["unit"],
        f'{second_var} low end': parametric_inputs[second_var]["low end"],
        f'{second_var} high end': parametric_inputs[second_var]["high end"],
        f'{second_var} step size': parametric_inputs[second_var]["step size"],
    }
    # create dictionary which will hold all the data points for this given parametric study
    param_results_dict = {
        f'{first_var}': [], # x-axis
        f'{second_var}': [], # y-axis
        'specific impulse': [], # z-axis
        'apogee reached': [], # a list of T/F values indicating whether the target apogee was reached with this input
        'rocket parameters for given input': [] # a list of dicts, containing details about the rocket for each iteration
    }
    # selectively remove rocket input variables
        # Mo and OF are mututally exclusive: if OF is selected, must remove Mo from rocket inputs
    remove_unnecessary_rocket_input_vars(rocket_input_values, first_var, second_var)
    
    # generate list of input values to traverse over
    x_vals = numpy.arange(param_iteration_dict[f'{first_var} low end'], param_iteration_dict[f'{first_var} high end'] + param_iteration_dict[f'{first_var} step size'], param_iteration_dict[f'{first_var} step size'])
    y_vals = numpy.arange(param_iteration_dict[f'{second_var} low end'], param_iteration_dict[f'{second_var} high end'] + param_iteration_dict[f'{second_var} step size'], param_iteration_dict[f'{second_var} step size'])
    
    # estimate runtime
    print(f"{blue}Performing parametric analysis on {first_var} and {second_var}: need to converge on target apogee {len(x_vals) * len(y_vals)} times{reset}") if simulation_settings["debug comments"] else None
    i=1
    # cartesian product traversal over the range of values given
    for x, y in itertools.product(x_vals, y_vals):
        # initialize rocket inputs dict
        rocket_input_values[first_var] = x
        rocket_input_values[second_var] = y
        # run fuel mass convergence algo
        fuel_convergence_results = fuel_mass_convergence.iterate_over_inner_radius(rocket_input_values, constants_dict, simulation_settings)
        # get Isp & rocket parameters
        rocket_parameters = fuel_convergence_results["rocket parameters"]
        Isp = rocket_parameters["Isp"]
        # add new values to param_results_dict
        param_results_dict[f'{first_var}'].append(x)
        param_results_dict[f'{second_var}'].append(y)
        param_results_dict["specific impulse"].append(Isp)
        param_results_dict["apogee reached"].append(rocket_parameters["apogee reached T/F"])
        param_results_dict["rocket parameters for given input"].append(rocket_parameters)
        
        # status report
        apogee_reached = rocket_parameters["apogee reached T/F"]
        color = green if apogee_reached else red
        print(f'{round(100 * i/(len(x_vals)*len(y_vals)),1)}% -- {blue}{first_var}: {x}, {second_var}: {y}, Isp: {Isp} -- reached apogee target: {color}{apogee_reached}{reset}') if simulation_settings["debug comments"] else None
        
        # save values and print status report if i is divisible by 100 (every 100 iterations)
        #if i % 100 == 0:
        #    print_param_hundredly_status_report(i, param_results_dict, round(100 * i/(len(x_vals)),1))
        #    output_handler.save_in_progress_param_simulation(param_iteration_dict, param_results_dict)
        
        i+=1
    
    return {
        "parametric study input dict": param_iteration_dict, 
        "parametric study output dict": param_results_dict}
    
# --------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY 3: Three Inputs --------------------
# --------------------------------------------------------------------------

def run_param_three_inputs(rocket_input_values, constants_dict, parametric_inputs, simulation_settings):
    # determine variable to parametrize
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    first_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]]
    second_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[1]]
    third_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[2]]
    # check that input variables are not mutually exclusive
    param_study_inputs_not_mutually_exclusive_checker(first_var, second_var, third_var)
    # create dictionary with bounds & step size
    param_iteration_dict = {
        "study name": parametric_inputs.get("study name"),
        f'{first_var} unit': parametric_inputs[first_var]["unit"],
        f'{first_var} low end': parametric_inputs[first_var]["low end"],
        f'{first_var} high end': parametric_inputs[first_var]["high end"],
        f'{first_var} step size': parametric_inputs[first_var]["step size"],
        f'{second_var} unit': parametric_inputs[second_var]["unit"],
        f'{second_var} low end': parametric_inputs[second_var]["low end"],
        f'{second_var} high end': parametric_inputs[second_var]["high end"],
        f'{second_var} step size': parametric_inputs[second_var]["step size"],
        f'{third_var} unit': parametric_inputs[third_var]["unit"],
        f'{third_var} low end': parametric_inputs[third_var]["low end"],
        f'{third_var} high end': parametric_inputs[third_var]["high end"],
        f'{third_var} step size': parametric_inputs[third_var]["step size"],
    }
    # create dictionary which will hold all the data points for this given parametric study
    param_results_dict = {
        f'{first_var}': [], # x-axis
        f'{second_var}': [], # y-axis
        f'{third_var}': [], # s-axis (will be a slider in the graph)
        'specific impulse': [], # z-axis
        'apogee reached': [], # a list of T/F values indicating whether the target apogee was reached with this input
        'rocket parameters for given input': [] # a list of dicts, containing details about the rocket for each iteration
    }
    # selectively remove rocket input variables
        # Mo and OF are mututally exclusive: if OF is selected, must remove Mo from rocket inputs
    remove_unnecessary_rocket_input_vars(rocket_input_values, first_var, second_var, third_var)
    
    # generate list of input values to traverse over
    x_vals = numpy.arange(param_iteration_dict[f'{first_var} low end'], param_iteration_dict[f'{first_var} high end'] + param_iteration_dict[f'{first_var} step size'], param_iteration_dict[f'{first_var} step size'])
    y_vals = numpy.arange(param_iteration_dict[f'{second_var} low end'], param_iteration_dict[f'{second_var} high end'] + param_iteration_dict[f'{second_var} step size'], param_iteration_dict[f'{second_var} step size'])
    s_vals = numpy.arange(param_iteration_dict[f'{third_var} low end'], param_iteration_dict[f'{third_var} high end'] + param_iteration_dict[f'{third_var} step size'], param_iteration_dict[f'{third_var} step size'])

    # estimate runtime
    print(f"{blue}Performing parametric analysis on {first_var}, {second_var}, and {third_var}: need to converge on target apogee {len(x_vals) * len(y_vals) * len(s_vals)} times{reset}") if simulation_settings["debug comments"] else None
    i=1 # keeps track iteration looping
    # cartesian product traversal over the range of values given
    for x, y, s in itertools.product(x_vals, y_vals, s_vals):
        # initialize rocket inputs dict
        rocket_input_values[first_var] = x
        rocket_input_values[second_var] = y
        rocket_input_values[third_var] = s
        # run fuel mass convergence algo
        fuel_convergence_results = fuel_mass_convergence.iterate_over_inner_radius(rocket_input_values, constants_dict, simulation_settings)
        # get Isp & rocket parameters
        rocket_parameters = fuel_convergence_results["rocket parameters"]
        Isp = rocket_parameters["Isp"]
        # add new values to param_results_dict
        param_results_dict[f'{first_var}'].append(x)
        param_results_dict[f'{second_var}'].append(y)
        param_results_dict[f'{third_var}'].append(s)
        param_results_dict["specific impulse"].append(Isp)
        param_results_dict["apogee reached"].append(rocket_parameters["apogee reached T/F"])
        param_results_dict["rocket parameters for given input"].append(rocket_parameters)
        
        # status report
        apogee_reached = rocket_parameters["apogee reached T/F"]
        color = green if apogee_reached else red
        print(f'{round(100 * i/(len(x_vals)*len(y_vals)*len(s_vals)),1)}% -- {blue}{first_var}: {x}, {second_var}: {y}, {third_var}: {s}, Isp: {Isp} -- reached apogee target: {color}{apogee_reached}{reset}') if simulation_settings["debug comments"] else None
        
        # save values and print status report if i is divisible by 100 (every 100 iterations)
        #if i % 100 == 0:
        #    print_param_hundredly_status_report(i, param_results_dict, round(100 * i/(len(x_vals)),1))
        #    output_handler.save_in_progress_param_simulation(param_iteration_dict, param_results_dict, parametric_inputs)
        
        i+=1
    
    return {
        "parametric study input dict": param_iteration_dict, 
        "parametric study output dict": param_results_dict}


# ---------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY HELPER FUNCTIONS --------------------
# ---------------------------------------------------------------------------

def param_study_inputs_not_mutually_exclusive_checker(first_var, second_var, third_var=None):
    # mutually exclusive pairs, contained in a list of sets
    mutually_exclusive_pairs = [
        {"oxidizer mass flow rate", "average oxidizer to fuel ratio"},
        {"var1", "var2"}, # expand list here
    ]
    # check if any inputs are identical (should never happen)
    inputs = [first_var, second_var]
    if third_var is not None:
        inputs.append(third_var)
    if len(inputs) != len(set(inputs)):
        raise ValueError(f"Inputs are not mutually exclusive: {inputs} contain identical values.")
    # check mutually exclusive pairs
    input_set = set(inputs)
    for pair in mutually_exclusive_pairs:
        if pair.issubset(input_set):
            raise ValueError(
                f"Inputs are not mutually exclusive: {pair} cannot be selected simultaneously."
            )
            
def remove_unnecessary_rocket_input_vars(rocket_input_values, first_var, second_var=None, third_var=None):
    # Mo and OF are mututally exclusive: if OF is selected, must remove Mo from rocket inputs
    if "average oxidizer to fuel ratio" in (first_var, second_var, third_var):
        #del rocket_input_values["oxidizer mass flow rate"]
        rocket_input_values.pop("oxidizer mass flow rate", None)

# every hundred parametric study iterations, print a report of how it's all going
def print_param_hundredly_status_report(nth_iteration, param_results_dict, percent_done):
    reached_apogee_TF_list = param_results_dict["apogee reached"]
    percentage_reached_apogee = 100 * sum(reached_apogee_TF_list)/len(reached_apogee_TF_list)
    print()
    print(f'{blue}Iteration {nth_iteration}: {percent_done}% of simulation finished. {percentage_reached_apogee}% of rockets have reached target apogee{reset}')
    print()