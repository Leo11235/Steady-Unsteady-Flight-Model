import steady_backend.graphing as graphing
import steady_backend.unit_conversion as convert
import numpy as np
from datetime import datetime
import os, json

from steady_backend.initialize_variables import shortened_variable_names as shortened_names

 


def print_simulation_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict):
    if simulation_type.get("optimize values for unsteady"):
        None
    elif simulation_type.get("simple fuel mass convergence"):
        # print fuel mass convergence results
        print_fuel_mass_convergence(rocket_inputs, simulation_type, simulation_settings, output_dict)
    elif simulation_type.get("parametric study"):
        print_parametric_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict)
    
        
        
# print results from simple fuel mass convergence
def print_fuel_mass_convergence(rocket_inputs, simulation_type, simulation_settings, output_dict):
    print_dict(output_dict)
    if simulation_settings["output unit system"] == "MRT":
        convert.convert_rocket_inputs_to_MRT(rocket_inputs) 
        convert.convert_rocket_parameters_to_MRT(output_dict["rocket parameters"])
    constants_dict = output_dict["constants dict"] if "constants dict" in output_dict else None
    rocket_parameters = output_dict["rocket parameters"] if "rocket parameters" in output_dict else None
    flight_data = output_dict["flight data"] if "flight data" in output_dict else None
    print("Rocket Inputs: ")
    print_dict(rocket_inputs)
    print("Rocket Parameters: ")
    print_dict(rocket_parameters)
    graphing.graph_everything(flight_data) if simulation_settings["show graphs"] else None
        
# print parametric study results
def print_parametric_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict):
    # print("Rocket inputs:")
    # print_dict(rocket_inputs)
    # print("Simulation type: ")
    # print_dict(simulation_type)
    # print("Simulation settings: ")
    # print_dict(simulation_settings)
    # print("Parametric inputs: ")
    # print_dict(parametric_inputs)
    # print("Outputs dict: ")
    # print_dict(output_dict)
    
    parametric_inputs_dict = output_dict["parametric study input dict"]
    parametric_outputs_dict = output_dict["parametric study output dict"]
    
    # convert np.float64 values in lists to normal floats (json doesn't support np.float64)
    if "parametric study output dict" in output_dict:
        output_dict["parametric study output dict"] = convert_numpy_types(output_dict["parametric study output dict"])
    
    # convert values to IMP as appropriate
    if simulation_settings["output unit system"] == "MRT": 
        convert.convert_rocket_inputs_to_MRT(rocket_inputs)
        convert.convert_parametric_inputs_to_MRT(parametric_inputs)
        convert.convert_param_output_dict_to_MRT(output_dict)
    
    # create a dictionary for graphing and for saving data, depending on simulation settings
    param_master_dict = {
            "title": generate_param_study_name(parametric_inputs),
            "parametric inputs": output_dict["parametric study input dict"],
            "parametric outputs": output_dict["parametric study output dict"]
        }
        
    graphing.graph_param(param_master_dict) if simulation_settings.get("show graphs") else None
    save_parametric_data(param_master_dict, rocket_inputs, parametric_inputs, parametric_inputs_dict, parametric_outputs_dict, simulation_settings) if simulation_settings.get("save simulation data") else None
    
    



# HELPER FUNCTIONS

# print a dictionary in a nice way
def print_dict(input_dictionary):
    print()
    max_key_length = max(len(str(key)) for key in input_dictionary.keys())
    for key, value in input_dictionary.items():
        fancyshmancy = "." * (max_key_length - len(str(key)) + 2) 
        print(f"{key} {fancyshmancy} {value}")
    print()

# take given rocket parameters, format them in the desired way to be viewed properly in the resulting graph
def format_parametric_study_rocket_parameters(rocket_parameters, rocket_inputs, simulation_settings):
    return {
        "Rocket": f'{rocket_inputs["rocket name"]}',
        "Specific Impulse": f'{round(rocket_parameters["Isp"], 2)} [s]',
        "Total Impulse": f'{round(rocket_parameters["total impulse"], 1)} [Ns]',
        "Burntime": f'{round(rocket_parameters["burntime"], 2)} [s]',
        
        "Thrust": f'{round(rocket_parameters["thrust"]/1000, 5)} [kN]',
        "Pad TtW": f'{round(rocket_parameters["thrust to weight ratio"], 3)}',
        
        "Wet Mass": f'{round(rocket_parameters["wet mass"], 2)} [kg]',
        "Ox Mass": f'{round(rocket_parameters["oxidizer mass"], 2)} [kg]',
        "Fuel Mass": f'{round(rocket_parameters["fuel mass"], 2)} [kg]',
        "Fuel Cell Inner Radius": f'{round(rocket_parameters["initial internal fuel radius"], 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}',
        
        "OF Ratio": f'{round(rocket_parameters["average oxidizer to fuel ratio"], 2)}' if "average oxidizer to fuel ratio" in rocket_parameters else None,
        "Chamber Temp": f'{round(rocket_parameters["chamber temperature"], 2)} [K]',
        
        "Nozzle Throat Area": f'{round(rocket_parameters["nozzle throat area"], 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle throat area" in rocket_parameters else None,
        "Nozzle Throat Radius": f'{round(rocket_parameters["nozzle throat radius"], 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle throat radius" in rocket_parameters else None,
        "Nozzle Outlet Area": f'{round(rocket_parameters["nozzle exit area"], 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle exit area" in rocket_parameters else None,
        "Nozzle Outlet Radius": f'{round(rocket_parameters["nozzle exit radius"], 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle exit radius" in rocket_parameters else None,
        "Nozzle Gas Exit Velocity": f'{round(rocket_parameters["nozzle gas exit velocity"], 2)} [m/s]',
        "Nozzle Gas Exit Temperature": f'{round(rocket_parameters["nozzle gas exit temperature"], 2)} [K]',
    }
   
# helper function to recursively convert numpy types to native Python types
def convert_numpy_types(obj):
    if isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(i) for i in obj]
    elif isinstance(obj, np.ndarray):  # If it's a numpy array, convert it to a list
        return obj.tolist()
    elif isinstance(obj, (np.float64, np.int64)):  # If it's a numpy scalar, convert to float or int
        return obj.item()
    else:
        return obj  # Return the object as is if it's not a numpy type
    

# save parametric study data
def save_parametric_data(param_master_dict, rocket_inputs, parametric_inputs, parametric_inputs_dict, parametric_outputs_dict, simulation_settings):
    # parametric data dict is a dict of dicts = 
        # "title": "title string generated in generate_study_name"
        # "parametric inputs": output_dict["parametric study input dict"],
        # "parametric outputs": output_dict["parametric study output dict"]
    
    # checks first if a file has already been generated for this parametric study, else, creates a new one
    # nvm not implemented yet
    
    # determine variables that were parametrized
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    first_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]] if len(indices_with_ones) > 0 else None
    second_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[1]] if len(indices_with_ones) > 1 else None
    third_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[2]] if len(indices_with_ones) > 2 else None
    
    variables = [var for var in [first_var, second_var, third_var] if var is not None]
    if len(variables) == 1:
        description_sentence = f"Optimized Isp as a function of {variables[0]} on range ({parametric_inputs_dict[f'{variables[0]} low end']}, {parametric_inputs_dict[f'{variables[0]} high end']} - step: {parametric_inputs_dict[f'{variables[0]} step size']})"
        description_range = f'({shortened_names[f'{variables[0]}']}[{parametric_inputs_dict[f'{variables[0]} low end']}, {parametric_inputs_dict[f'{variables[0]} high end']}])'
        gaph_title = f'{list(parametric_outputs_dict.keys())[1].capitalize()} vs {list(parametric_outputs_dict.keys())[0].capitalize()}'
    elif len(variables) == 2:
        description_sentence = f"Optimized Isp as a function of {variables[0]} and {variables[1]}"
        description_range = f'({shortened_names[f'{variables[0]}']}[{parametric_inputs_dict[f'{variables[0]} low end']}, {parametric_inputs_dict[f'{variables[0]} high end']}], {shortened_names[f'{variables[1]}']}[{parametric_inputs_dict[f'{variables[1]} low end']}, {parametric_inputs_dict[f'{variables[1]} high end']}])'
        gaph_title = f'{list(parametric_outputs_dict.keys())[2].capitalize()} vs {list(parametric_outputs_dict.keys())[0].capitalize()} and {list(parametric_outputs_dict.keys())[1].capitalize()}'
    elif len(variables) == 3:
        description_sentence = f"Optimized Isp as a function of {variables[0]}, {variables[1]}, and {variables[2]}"
        description_range = f'({shortened_names[f'{variables[0]}']}[{parametric_inputs_dict[f'{variables[0]} low end']}, {parametric_inputs_dict[f'{variables[0]} high end']}], {shortened_names[f'{variables[1]}']}[{parametric_inputs_dict[f'{variables[1]} low end']}, {parametric_inputs_dict[f'{variables[1]} high end']}], {shortened_names[f'{variables[2]}']}[{parametric_inputs_dict[f'{variables[2]} low end']}, {parametric_inputs_dict[f'{variables[2]} high end']}])'
        gaph_title = f'{list(parametric_outputs_dict.keys())[3].capitalize()} vs {list(parametric_outputs_dict.keys())[0].capitalize()}, {list(parametric_outputs_dict.keys())[1].capitalize()}, and {list(parametric_outputs_dict.keys())[2].capitalize()}'
    else:
        description_sentence = "Optimized Isp"
    
    now = datetime.now()
    date_str = now.strftime("%d-%m-%Y")
    time_str = now.strftime("%H-%M-%S") 
    
    # define output folder & make sure it exists
    folder_path = os.path.join("static", "saved_graphs")
    os.makedirs(folder_path, exist_ok=True)
    file_name = f'{rocket_inputs["rocket name"]} -- {description_range} -- ({date_str}, {time_str})'
    # Replace invalid characters in file name if needed
    file_name = file_name.replace(":", "-").replace("/", "-").replace("\\", "-")
    # Define the full file path
    file_path = os.path.join(folder_path, file_name+".json")
    
    # save to JSON
    try:
        with open(file_path, "w") as f:
            json.dump(param_master_dict, f, indent=4)
        print(f"Data saved to {file_path}") if simulation_settings["debug comments"] == "comments" or simulation_settings["debug comments"] != "none" else None
    except PermissionError:
        print(f"Permission denied: Unable to write to {file_path}")
    except Exception as e:
        print(f"An error occurred while saving data: {e}")
    
    print("File successfully saved!")

    
# returns a string which is the name of the parametric study
def generate_param_study_name(parametric_inputs):
    # determine variables that were parametrized
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    first_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]] if len(indices_with_ones) > 0 else None
    second_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[1]] if len(indices_with_ones) > 1 else None
    third_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[2]] if len(indices_with_ones) > 2 else None
    return f"Isp vs {first_var.capitalize() if first_var != None else ""}{(f', {second_var.capitalize()}') if second_var != None else ""}{(f', and {third_var.capitalize()}') if third_var != None else ""}"