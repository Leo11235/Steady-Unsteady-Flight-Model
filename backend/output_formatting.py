# convert units for frontend use
def format_output(rocket_inputs, rocket_parameters):
    rocket_inputs["chamber pressure"] = rocket_inputs["chamber pressure"] * 0.000145038 # Pa to psi
    rocket_inputs["fuel external diameter"] = rocket_inputs["fuel external diameter"] * 39.3701 # m to in
    rocket_inputs["fuel length"] = rocket_inputs["fuel length"] * 39.3701 # m to in
    rocket_inputs["target apogee"] = rocket_inputs["target apogee"] * 3.28084 # m to ft
    rocket_inputs["launch site altitude"] = rocket_inputs["launch site altitude"] * 3.28084 # m to ft
    rocket_inputs["rocket external diameter"] = rocket_inputs["rocket external diameter"] * 39.3701 # m to in
    
    rocket_parameters["initial internal fuel radius"] = rocket_parameters["initial internal fuel radius"] * 39.3701 # m to in
    rocket_parameters["nozzle throat area"] = rocket_parameters["nozzle throat area"] * 1550 # m^2 to in^2
    rocket_parameters["nozzle throat radius"] = rocket_parameters["nozzle throat radius"] * 39.3701 # m to in
    rocket_parameters["nozzle gas exit pressure"] = rocket_parameters["nozzle gas exit pressure"] * 0.000145038 # Pa to psi
    rocket_parameters["nozzle exit area"] = rocket_parameters["nozzle exit area"] * 1550 # m^2 to in^2
    rocket_parameters["nozzle exit radius"] = rocket_parameters["nozzle exit radius"] * 39.3701 # m to in

    return rocket_inputs, rocket_parameters
