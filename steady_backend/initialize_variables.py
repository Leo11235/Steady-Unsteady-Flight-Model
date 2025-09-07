import json
import steady_backend.unit_conversion as convert_units


# dictionary containing universal constants used throughout the program
constants_dict = {
    # universal constants
    "universal gas constant": 8.31446261815324, # J / (mol * K)
    
    # planet constants (Earth comes as a preset, can be written over)
    "atmosphere": True, # body has an atmosphere
    "sea level gravity": 9.80665, # m/s^2
    "ambient sea level atmospheric pressure": 101325,  # Pa 
    "planet radius": 6378000, # m
    "sea level air density": 1.225, # kg/m^3
    "sea level temperature": 288.15, # K
    "stratosphere temperature": 216.65, # K
    "air molar mass": 0.0289644, # kg/mol
    "temperature lapse rate in the troposphere": 0.0065, # K/m

    # for computation purposes
    "number of timesteps in burn": 100, # for calculating rocket ascent
    "tolerated apogee difference": .1, # m, allowed space between target and real apogee
    "smallest allowed inner fuel radius": 0.01, # m; equivalent to 1 cm, 0.0127 suggested as well (0.5 in)
}

# dict containing all the kinds of errors that may occur in the program
# this gets passed to the frontend so the user can see why their inputs failed
# errors must be formatted in this way: 
# "errortype error": ...
# "errortype error information": ...
error_dict = {
    "error": False, # if this is still set to false when it gets passed to the frontend, this entire dictionary is ignored
    "bisection error": False, # if the rocket is unable to reach target apogee
    "bisection error information": "",
    "unknown error": False, 
    "unknown error information": "",
}

shortened_variable_names = { # shortened names of variables you can run parametric studies on, this is used for file creation
        "oxidizer mass flow rate": "Mo",
        "average oxidizer to fuel ratio": "OF",
        "chamber pressure": "CP",
        "fuel length": "Lf",
        "target apogee": "Target Apogee",
        "launch angle": "Launch Angle",
        "drag coefficient": "Drag Coefficient", 
        "dry mass": "Dry Mass",
        "rocket external diameter": "Re",
        "fuel grain density": "pf",
        "nozzle throat radius": "Rtn",
        "nozzle exit radius": "Ren",
    }

variable_units_SI = { 
    "oxidizer mass flow rate": "kg/s",
    "average oxidizer to fuel ratio": "dimensionless",
    "chamber pressure": "Pa",
    "nozzle throat radius": "m",
    "nozzle exit radius": "m",
    "fuel length": "m",    
}

def initialize_constants_dict(planet):
    # load JSON containing planet info
    planets_json = "./static/planets_data/Extended_Planet_Attributes.json"
    with open(planets_json, 'r') as json_file:
        planet_attributes = json.load(json_file)
        
    planet_dict = planet_attributes[planet] # replace the planet string (ie "Earth") with a dict containing its attributes       
        
    # overwrite planet attributes
    constants_dict["atmosphere"] = planet_dict["atmosphere"]
    constants_dict["sea level gravity"] = planet_dict["sea level gravity"]
    constants_dict["planet radius"] = planet_dict["planet radius"]
    if constants_dict["atmosphere"]:
        constants_dict["sea level air density"] = planet_dict["sea level air density"]
        constants_dict["sea level temperature"] = planet_dict["sea level temperature"]
        constants_dict["air molar mass"] = planet_dict["air molar mass"]
    if planet == "Earth": # just to make air density calculations extra precise
        constants_dict["ambient sea level atmospheric pressure"] = planet_dict["ambient sea level atmospheric pressure"]
        constants_dict["stratosphere temperature"] = planet_dict["stratosphere temperature"]
        constants_dict["troposphere temperature lapse rate"] = planet_dict["troposphere temperature lapse rate"]
    
    return constants_dict

# create a rocket_input_values dict containing information about the rocket
def initialize_rocket_input_values_dict(rocket_inputs):
    convert_units.convert_to_SI(rocket_inputs)
    
    rocket_input_values = {
        "rocket name": rocket_inputs["rocket name"],
        "planet": rocket_inputs["planet"],
        "target apogee": rocket_inputs["target apogee"][0],
        "oxidizer mass flow rate": rocket_inputs["oxidizer mass flow rate"][0],
        "chamber pressure": rocket_inputs["chamber pressure"][0], 
        "fuel external diameter": rocket_inputs["fuel external diameter"][0],
        "fuel length": rocket_inputs["fuel length"][0],
        "launch site altitude": rocket_inputs["launch site altitude"][0],
        "dry mass": rocket_inputs["dry mass"][0],
        "rocket external diameter": rocket_inputs["rocket external diameter"][0],
        "fuel grain density": rocket_inputs["fuel grain density"][0],
        "liquid oxidizer type": rocket_inputs["liquid oxidizer type"],
        "solid fuel type": rocket_inputs["solid fuel type"],
        "launch angle": rocket_inputs["launch angle"],
        "drag coefficient": rocket_inputs["drag coefficient"],
        "regression rate exponent": rocket_inputs["regression rate exponent"],
        "regression rate scaling coefficient": rocket_inputs["regression rate scaling coefficient"],
    }
    
    return rocket_input_values