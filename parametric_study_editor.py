# change the bounds and / or step size of an existing parametric study
# takes as input an existing parametric study json file and a dict containing the details to be edited, returns a new parametric study json
# can be used to expand or crop the bounds of a study, give an existing study a higher 'resolution', etc.

import json
import numpy as np
import steady_backend.unit_conversion as convert

file_name = "Ancalagon 2 -- (CP[100, 1000], Lf[7.87, 78.74]) -- (30-01-2025, 15-01-35).json"
file_path = f"./static/saved_graphs/{file_name}" + (".json" if not file_name.endswith(".json") else "")

edits_dict = { # contains new range and stepsize to parametrize over, any fields not included will default to what already exists; also contains instructions of what to do with new results
            "variables to redefine": ["chamber pressure"],
            "chamber pressure unit": "IMP",
            "chamber pressure low end": 400,
            "chamber pressure high end": 800,
            "chamber pressure step size": 50,
            
            "graph results": True,
            "save results": True,
}

def edit_parametric_study(P_Study_json, edits_dict):
    # unpack json
    with open(P_Study_json, 'r') as json_file:
        parametric_data = json.load(json_file)
    PS_description = parametric_data["description"]
    PS_title = parametric_data["title"]
    PS_inputs = parametric_data["parametric inputs"]
    PS_outputs = parametric_data["parametric outputs"]
    
    # convert any IMP to SI values
    print(edits_dict)
    convert.convert_param_to_SI(edits_dict)
    print(edits_dict)

    # print(PS_description)
    # print(PS_title)
    print(PS_inputs)
    # print(PS_outputs)
    
    # define ranges to parametrize over
    vars_to_redefine = 0 # number of variables getting parametrized over
    vars_to_redefine_dict = {} # dictionary containing lists of values to parametrize over
    for var in edits_dict["variables to redefine"]:
        vars_to_redefine_dict[var] = PS_inputs[f'{var} unit'] # returns output values in the same measurement system as input
        # for each variable, create a list of the values to parametrize over: 
        low_end = PS_inputs[f'{var} low end']
        vars_to_redefine_dict[f'{var} list'] = np.arange(low_end, PS_inputs[f'{var} high end'], PS_inputs[f'{var} step size'])
        vars_to_redefine += 1
    
    print(vars_to_redefine_dict)
    
    # remove any values already inclued (+- 0.1% to account for conversion errors)



if __name__ == "__main__":
    edit_parametric_study(file_path, edits_dict)
    
    