import steady_backend.graphing as graphing
import steady_backend.unit_conversion as convert
from steady_backend.initialize_variables import shortened_variable_names as shortened_names
import matplotlib.pyplot as plt
import os, json
from datetime import datetime
import numpy as np
    
    
def print_simulation_results(rocket_inputs, simulation_type, simulation_settings, parametric_inputs, output_dict):
    # convert np.float64 values in lists to normal floats
    if "parametric study output dict" in output_dict:
        output_dict["parametric study output dict"] = convert_numpy_types(output_dict["parametric study output dict"])
    # convert values to the desired output unit system
    if simulation_settings["output unit system"] == "MRT":
        convert.convert_rocket_inputs_to_MRT(rocket_inputs)
        convert.convert_parametric_inputs_to_MRT(parametric_inputs)
        convert.convert_output_dict_to_MRT(output_dict)
    # print output according to simulation type & settings
    if simulation_type["optimize values for unsteady"]:
        None # not yet implemented
    
    elif simulation_type["simple fuel mass convergence"]:
        rocket_inputs = output_dict["rocket inputs"] if "rocket inputs" in output_dict else None
        constants_dict = output_dict["constants dict"] if "constants dict" in output_dict else None
        rocket_parameters = output_dict["rocket parameters"] if "rocket parameters" in output_dict else None
        flight_data = output_dict["flight data"] if "flight data" in output_dict else None
        print_dict(rocket_inputs)
        print_dict(constants_dict)
        print_dict(rocket_parameters)
        graphing.graph_everything(flight_data) if simulation_settings["show graphs"] else None
        
    elif simulation_type["parametric study"]: # and sum(parametric_inputs["variables to parametrize"]) == 1: 
        parametric_inputs_dict = output_dict["parametric study input dict"]
        parametric_outputs_dict = output_dict["parametric study output dict"]
        for index, item in enumerate(parametric_outputs_dict["rocket parameters for given input"]):
            parametric_outputs_dict["rocket parameters for given input"][index] = format_parametric_study_rocket_parameters(item, rocket_inputs, simulation_settings)
        # save and graph data, if applicable for each
        save_parametric_data(rocket_inputs, parametric_inputs, parametric_inputs_dict, parametric_outputs_dict, simulation_settings) if simulation_settings["save simulation data"] else None
        graphing.graph_param({"title": f'{list(parametric_outputs_dict.keys())[1]} vs {list(parametric_outputs_dict.keys())[0]}', 
                                        "parametric inputs": parametric_inputs_dict,
                                        "parametric outputs": parametric_outputs_dict}) if simulation_settings["show graphs"] else None

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
        "Fuel Cell Inner Radius": f'{round(rocket_parameters["initial internal fuel radius"] * (39.3701 if simulation_settings["output unit system"] == "MRT" else 100), 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}',
        # fuel, ox, fuel+ox mass flow rates?
        
        "OF Ratio": f'{round(rocket_parameters["average oxidizer to fuel ratio"], 2)}' if "average oxidizer to fuel ratio" in rocket_parameters else None,
        "Chamber Temp": f'{round(rocket_parameters["chamber temperature"], 2)} [K]',
        
        "Nozzle Throat Area": f'{round(rocket_parameters["nozzle throat area"] * (1550 if simulation_settings["output unit system"] == "MRT" else 10000), 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle throat area" in rocket_parameters else None,
        "Nozzle Throat Radius": f'{round(rocket_parameters["nozzle throat radius"] * (39.3701 if simulation_settings["output unit system"] == "MRT" else 100), 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle throat radius" in rocket_parameters else None,
        "Nozzle Outlet Area": f'{round(rocket_parameters["nozzle exit area"] * (1550 if simulation_settings["output unit system"] == "MRT" else 10000), 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle exit area" in rocket_parameters else None,
        "Nozzle Outlet Radius": f'{round(rocket_parameters["nozzle exit radius"] * (39.3701 if simulation_settings["output unit system"] == "MRT" else 100), 4)} {"[in]" if simulation_settings["output unit system"] == "MRT" else "[cm]"}' if "nozzle exit radius" in rocket_parameters else None,
        "Nozzle Gas Exit Velocity": f'{round(rocket_parameters["nozzle gas exit velocity"], 2)} [m/s]',
        "Nozzle Gas Exit Temperature": f'{round(rocket_parameters["nozzle gas exit temperature"], 2)} [K]',
    }
    
# save parametric data in a nely created json file
def save_parametric_data(rocket_inputs, parametric_inputs, parametric_inputs_dict, parametric_outputs_dict, simulation_settings):
    # determine variables that were parametrized
    indices_with_ones = [i for i, value in enumerate(parametric_inputs["variables to parametrize"]) if value == 1] # list of indices with values set to 1 (meaning the corresponding variable is to be parametrized)
    first_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[0]] if len(indices_with_ones) > 0 else None
    second_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[1]] if len(indices_with_ones) > 1 else None
    third_var = parametric_inputs["options of variables to parametrize"][indices_with_ones[2]] if len(indices_with_ones) > 2 else None
    
    # dynamically create file name and description, and get the date and time
    # also dynamically create abbreviated names for variables and generate a kinda statement to describe their ranges
    
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
    
    # create json data
    data = {
        "description": 
            f'{rocket_inputs["rocket name"]} launched from {rocket_inputs["planet"]} with target apogee {round(rocket_inputs["target apogee"][0],1)} meters ({round(rocket_inputs["target apogee"][0]*3.28084,1)} ft). {description_sentence}. '
            f'Created on {date_str} at {time_str}',
        "title": gaph_title,
        "parametric inputs": parametric_inputs_dict,
        "parametric outputs": parametric_outputs_dict
    }
    
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
            json.dump(data, f, indent=4)
        print(f"Data saved to {file_path}") if simulation_settings["debug comments"] == "comments" or simulation_settings["debug comments"] != "none" else None
    except PermissionError:
        print(f"Permission denied: Unable to write to {file_path}")
    except Exception as e:
        print(f"An error occurred while saving data: {e}")
    
    
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