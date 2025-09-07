# convert any units from IMP to SI
def convert_to_SI(input_dict): # used for rocket inputs
    # rocket inputs
    if "target apogee" in input_dict and input_dict["target apogee"][1] == "IMP":
        input_dict["target apogee"][0] = float(input_dict["target apogee"][0]) * 0.3048 # ft to m
        input_dict["target apogee"][1] = "SI"
    if "oxidizer mass flow rate" in input_dict and input_dict["oxidizer mass flow rate"][1] == "IMP":
        input_dict["oxidizer mass flow rate"][0] = float(input_dict["oxidizer mass flow rate"][0]) * 0.453592 # lb/s to kg/s
        input_dict["oxidizer mass flow rate"][1] = "SI"
    if "chamber pressure" in input_dict and input_dict["chamber pressure"][1] == "IMP":
        input_dict["chamber pressure"][0] = float(input_dict["chamber pressure"][0]) * 6894.76 # psi to Pa
        input_dict["chamber pressure"][1] = "SI"
    if "fuel external diameter" in input_dict and input_dict["fuel external diameter"][1] == "IMP":
        input_dict["fuel external diameter"][0] = float(input_dict["fuel external diameter"][0]) * 0.0254 # in to m
        input_dict["fuel external diameter"][1] = "SI"
    if "fuel length" in input_dict and input_dict["fuel length"][1] == "IMP":
        input_dict["fuel length"][0] = float(input_dict["fuel length"][0]) * 0.0254 # in to m
        input_dict["fuel length"][1] = "SI"
    if "launch site altitude" in input_dict and input_dict["launch site altitude"][1] == "IMP":
        input_dict["launch site altitude"][0] = float(input_dict["launch site altitude"][0]) * 0.3048 # ft to m
        input_dict["launch site altitude"][1] = "SI"
    if "dry mass" in input_dict and input_dict["dry mass"][1] == "IMP":
        input_dict["dry mass"][0] = float(input_dict["dry mass"][0]) * 0.453592 # lb to kg
        input_dict["dry mass"][1] = "SI"
    if "rocket external diameter" in input_dict and input_dict["rocket external diameter"][1] == "IMP":
        input_dict["rocket external diameter"][0] = float(input_dict["rocket external diameter"][0]) * 0.0254 # in to m
        input_dict["rocket external diameter"][1] = "SI"
    if "fuel grain density" in input_dict and input_dict["fuel grain density"][1] == "IMP":
        input_dict["fuel grain density"][0] = float(input_dict["fuel grain density"][0]) * 16.0185 # lb/ft^3 to kg/m^3
        input_dict["fuel grain density"][1] = "SI"
    # ensure remaining values are floats
    input_dict["launch angle"] = float(input_dict["launch angle"]) if "launch angle" in input_dict else None
    input_dict["drag coefficient"] = float(input_dict["drag coefficient"]) if "drag coefficient" in input_dict else None
    input_dict["regression rate exponent"] = float(input_dict["regression rate exponent"]) if "regression rate exponent" in input_dict else None
    input_dict["regression rate scaling coefficient"] = float(input_dict["regression rate scaling coefficient"]) if "regression rate scaling coefficient" in input_dict else None
    
    '''
    if "unit" in input_dict and input_dict["unit"][1] == "IMP":
        input_dict["unit"][0] = float(input_dict["unit"][0]) ########
        input_dict["unit"][1] = "SI"
    
    if "value" in input_dict:
        input_dict["unit"] = float(input_dict["unit"])
    '''    
    return input_dict
    
# convert parametric study input values to SI
def convert_param_to_SI(param_dict):
    print(param_dict)
    if "oxidizer mass flow rate" in param_dict and param_dict["oxidizer mass flow rate"]["unit"] == "IMP":
        param_dict["oxidizer mass flow rate"]["low end"] = float(param_dict["oxidizer mass flow rate"]["low end"]) * 0.453592 # lb/s to kg/s
        param_dict["oxidizer mass flow rate"]["high end"] = float(param_dict["oxidizer mass flow rate"]["high end"]) * 0.453592 # lb/s to kg/s
        param_dict["oxidizer mass flow rate"]["step size"] = float(param_dict["oxidizer mass flow rate"]["step size"]) * 0.453592 # lb/s to kg/s
    if "chamber pressure" in param_dict and param_dict["chamber pressure"]["unit"] == "IMP":
        param_dict["chamber pressure"]["low end"] = float(param_dict["chamber pressure"]["low end"]) * 6894.76 # psi to Pa
        param_dict["chamber pressure"]["high end"] = float(param_dict["chamber pressure"]["high end"]) * 6894.76 # psi to Pa
        param_dict["chamber pressure"]["step size"] = float(param_dict["chamber pressure"]["step size"]) * 6894.76 # psi to Pa
    if "fuel length" in param_dict and param_dict["fuel length"]["unit"] == "IMP":
        param_dict["fuel length"]["low end"] = float(param_dict["fuel length"]["low end"]) * 0.0254 # in to m
        param_dict["fuel length"]["high end"] = float(param_dict["fuel length"]["high end"]) * 0.0254 # in to m
        param_dict["fuel length"]["step size"] = float(param_dict["fuel length"]["step size"]) * 0.0254 # in to m
    if "target apogee" in param_dict and param_dict["target apogee"]["unit"] == "IMP":
        param_dict["target apogee"]["low end"] = float(param_dict["target apogee"]["low end"]) * 0.3048 # ft to m
        param_dict["target apogee"]["high end"] = float(param_dict["target apogee"]["high end"]) * 0.3048 # ft to m
        param_dict["target apogee"]["step size"] = float(param_dict["target apogee"]["step size"]) * 0.3048 # ft to m
    if "dry mass" in param_dict and param_dict["dry mass"]["unit"] == "IMP":
        param_dict["dry mass"]["low end"] = float(param_dict["dry mass"]["low end"]) * 0.453592 # lb to kg
        param_dict["dry mass"]["high end"] = float(param_dict["dry mass"]["high end"]) * 0.453592 # lb to kg
        param_dict["dry mass"]["step size"] = float(param_dict["dry mass"]["step size"]) * 0.453592 # lb to kg
    if "rocket external diameter" in param_dict and param_dict["rocket external diameter"]["unit"] == "IMP":
        param_dict["rocket external diameter"]["low end"] = float(param_dict["rocket external diameter"]["low end"]) * 0.0254 # in to m
        param_dict["rocket external diameter"]["high end"] = float(param_dict["rocket external diameter"]["high end"]) * 0.0254 # in to m
        param_dict["rocket external diameter"]["step size"] = float(param_dict["rocket external diameter"]["step size"]) * 0.0254 # in to m
    if "fuel grain density" in param_dict and param_dict["fuel grain density"]["unit"] == "IMP":
        param_dict["fuel grain density"]["low end"] = float(param_dict["fuel grain density"]["low end"]) * 16.0185 # lb/ft^3 to kg/m^3
        param_dict["fuel grain density"]["high end"] = float(param_dict["fuel grain density"]["high end"]) * 16.0185 # lb/ft^3 to kg/m^3
        param_dict["fuel grain density"]["step size"] = float(param_dict["fuel grain density"]["step size"]) * 16.0185 # lb/ft^3 to kg/m^3
    
    
    # not implemented into the rest of the program yet, still good to have here
    if "nozzle throat radius" in param_dict and param_dict["nozzle throat radius"]["unit"] == "IMP":
        param_dict["nozzle throat radius"]["low end"] = float(param_dict["nozzle throat radius"]["low end"]) * 0.0254 # in to m
        param_dict["nozzle throat radius"]["high end"] = float(param_dict["nozzle throat radius"]["high end"]) * 0.0254 # in to m
        param_dict["nozzle throat radius"]["step size"] = float(param_dict["nozzle throat radius"]["step size"]) * 0.0254 # in to m
    if "nozzle exit radius" in param_dict and param_dict["nozzle exit radius"]["unit"] == "IMP":
        param_dict["nozzle exit radius"]["low end"] = float(param_dict["nozzle exit radius"]["low end"]) * 0.0254 # in to m
        param_dict["nozzle exit radius"]["high end"] = float(param_dict["nozzle exit radius"]["high end"]) * 0.0254 # in to m
        param_dict["nozzle exit radius"]["step size"] = float(param_dict["nozzle exit radius"]["step size"]) * 0.0254 # in to m




# -----------------------------------------
# ----------- OUTPUT CONVERSION -----------
# -----------------------------------------

# convert any values to MRT SI-IMP mishmash
def convert_rocket_inputs_to_MRT(rocket_inputs):
    rocket_inputs["target apogee"][0] *= 3.28084 # m to ft
    rocket_inputs["target apogee"][1] = "IMP"
    rocket_inputs["chamber pressure"][0] *= 0.000145038 # Pa to psi
    rocket_inputs["chamber pressure"][1] = "IMP"
    rocket_inputs["fuel external diameter"][0] *= 39.3701 # m to in
    rocket_inputs["fuel external diameter"][1] = "IMP"
    rocket_inputs["fuel length"][0] *= 39.3701 # m to in
    rocket_inputs["fuel length"][1] = "IMP"
    rocket_inputs["launch site altitude"][0] *= 3.28084 # m to ft
    rocket_inputs["launch site altitude"][1] = "IMP"
    rocket_inputs["rocket external diameter"][0] *= 3.28084 # m to ft
    rocket_inputs["rocket external diameter"][1] = "IMP"

# convert output rocket parameters dict to weird MRT mishmash of SI and IMP values
def convert_rocket_parameters_to_MRT(rocket_parameters):
    if "initial internal fuel radius" in rocket_parameters:
        rocket_parameters["initial internal fuel radius"] *= 39.3701 # m to in
    if "nozzle throat area" in rocket_parameters:
        rocket_parameters["nozzle throat area"] *= 1550  # m^2 to in^2
    if "nozzle throat radius" in rocket_parameters:
        rocket_parameters["nozzle throat radius"] *= 39.3701  # m to in
    if "nozzle gas exit pressure" in rocket_parameters:
        rocket_parameters["nozzle gas exit pressure"] *= 0.000145038  # Pa to psi
    if "nozzle exit area" in rocket_parameters:
        rocket_parameters["nozzle exit area"] *= 1550  # m^2 to in^2
    if "nozzle exit radius" in rocket_parameters:
        rocket_parameters["nozzle exit radius"] *= 39.3701  # m to in
    if "reached apogee" in rocket_parameters:
        rocket_parameters["reached apogee"] *= 3.28084  # m to ft
    
def convert_parametric_inputs_to_MRT(parametric_inputs):
    if "chamber pressure" in parametric_inputs:
        parametric_inputs["chamber pressure"]["unit"] = "IMP"
        parametric_inputs["chamber pressure"]["low end"] *= 0.000145038 # Pa to psi
        parametric_inputs["chamber pressure"]["high end"] *= 0.000145038 # Pa to psi
        parametric_inputs["chamber pressure"]["step size"] *= 0.000145038 # Pa to psi
    if "nozzle throat radius" in parametric_inputs:
        parametric_inputs["nozzle throat radius"]["unit"] = "IMP"
        parametric_inputs["nozzle throat radius"]["low end"] *= 39.3701 # m to in
        parametric_inputs["nozzle throat radius"]["high end"] *= 39.3701 # m to in
        parametric_inputs["nozzle throat radius"]["step size"] *= 39.3701 # m to in
    if "" in parametric_inputs:
        parametric_inputs["nozzle exit radius"]["unit"] = "IMP"
        parametric_inputs["nozzle exit radius"]["low end"] *= 39.3701 # m to in
        parametric_inputs["nozzle exit radius"]["high end"] *= 39.3701 # m to in
        parametric_inputs["nozzle exit radius"]["step size"] *= 39.3701 # m to in
    if "fuel length" in parametric_inputs:
        parametric_inputs["fuel length"]["unit"] = "IMP"
        parametric_inputs["fuel length"]["low end"] *= 39.3701 # m to in
        parametric_inputs["fuel length"]["high end"] *= 39.3701 # m to in
        parametric_inputs["fuel length"]["step size"] *= 39.3701 # m to in
    if "target apogee" in parametric_inputs:
        parametric_inputs["target apogee"]["unit"] = "IMP"
        parametric_inputs["target apogee"]["low end"] *= 3.28084 # m to ft
        parametric_inputs["target apogee"]["high end"] *= 3.28084 # m to ft
        parametric_inputs["target apogee"]["step size"] *= 3.28084 # m to ft
    if "rocket external diameter" in parametric_inputs:
        parametric_inputs["rocket external diameter"]["unit"] = "IMP"
        parametric_inputs["rocket external diameter"]["low end"] *= 39.3701 # m to in
        parametric_inputs["rocket external diameter"]["high end"] *= 39.3701 # m to in
        parametric_inputs["rocket external diameter"]["step size"] *= 39.3701 # m to in
    
    
def convert_param_output_dict_to_MRT(output_dict):
    param_inputs = output_dict["parametric study input dict"]
    param_outputs = output_dict["parametric study output dict"]
    if "chamber pressure low end" in param_inputs:
        param_inputs["chamber pressure low end"] *= 0.000145038 # Pa to psi
        param_inputs["chamber pressure high end"] *= 0.000145038 # Pa to psi
        param_inputs["chamber pressure step size"] *= 0.000145038 # Pa to psi
    if "nozzle throat radius low end" in param_inputs:
        param_inputs["nozzle throat radius low end"] *= 39.3701 # m to in
        param_inputs["nozzle throat radius high end"] *= 39.3701 # m to in
        param_inputs["nozzle throat radius step size"] *= 39.3701 # m to in
    if "nozzle exit radius low end" in param_inputs:
        param_inputs["nozzle exit radius low end"] *= 39.3701 # m to in
        param_inputs["nozzle exit radius high end"] *= 39.3701 # m to in
        param_inputs["nozzle exit radius step size"] *= 39.3701 # m to in
    if "fuel length low end" in param_inputs:
        param_inputs["fuel length low end"] *= 39.3701 # m to in
        param_inputs["fuel length high end"] *= 39.3701 # m to in
        param_inputs["fuel length step size"] *= 39.3701 # m to in
    if "target apogee low end" in param_inputs:
        param_inputs["target apogee low end"] *= 3.28084 # m to ft
        param_inputs["target apogee high end"] *= 3.28084 # m to ft
        param_inputs["target apogee step size"] *= 3.28084 # m to ft
    if "rocket external diameter low end" in param_inputs:
        param_inputs["rocket external diameter low end"] *= 39.3701 # m to in
        param_inputs["rocket external diameter high end"] *= 39.3701 # m to in
        param_inputs["rocket external diameter step size"] *= 39.3701 # m to in
    
    if "chamber pressure" in param_outputs:
        param_outputs["chamber pressure"] = [item*0.000145038 for item in param_outputs["chamber pressure"]] # Pa to psi
    if "nozzle throat radius" in param_outputs:
        param_outputs["nozzle throat radius"] = [item*39.3701 for item in param_outputs["nozzle throat radius"]] # m to in
    if "nozzle exit radius" in param_outputs:
        param_outputs["nozzle exit radius"] = [item*39.3701 for item in param_outputs["nozzle exit radius"]] # m to in
    if "fuel length" in param_outputs:
        param_outputs["fuel length"] = [item*39.3701 for item in param_outputs["fuel length"]] # m to in
    if "target apogee" in param_outputs:
        param_outputs["target apogee"] = [item*3.28084 for item in param_outputs["target apogee"]] # m to ft
    if "rocket external diameter" in param_outputs:
        param_outputs["rocket external diameter"] = [item*39.3701 for item in param_outputs["fuel length"]] # m to in
    
    # convert list of dicts (each dict is rocket parameters, we call that function on each one here)
    for RP in param_outputs["rocket parameters for given input"]:
        convert_rocket_parameters_to_MRT(RP)
