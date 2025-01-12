'''
this class handless the input file and format is in a way
that is nice to use for the rest of the program
'''

import json

# list of required inputs
required_keys = [ # no longer checked / needed, keeping for now for testing
    "rocket name",

    # combustion stuff
    "oxidizer mass flow rate",
    "chamber pressure",
    "fuel external diameter",
    "fuel length",

    # rocket info
    "target apogee",
    "launch site altitude",
    "dry mass",
    "rocket external diameter",
    "drag coefficient",
    "launch angle",

    # idk yet
    "prop dry mass",
    "prop diameter",
    "prop length",

    # for PROPEP
    "liquid oxidizer type",
    "oxidizer molecule code",
    "solid fuel type",
    "fuel molecule code",

    # from research
    "fuel grain density",
    "regression rate scaling coefficient",
    "regression rate exponent",
]

# dictionary containing universal constants
# may end up adding lots of random stuff here
constants_dict = {
    # natural constants
    "sea level gravity": 9.80665, # m/s^2
    "universal gas constant": 8.31446261815324, # J / (mol * K)
    "ambient sea level atmospheric pressure": 101325,  # Pa 
    "earth radius": 6378000, # m
    "sea level air density": 1.225, # kg/m^3
    "sea level temperature": 288.15, # K
    "stratosphere temperature": 216.65, # K
    "air molar mass": 0.0289644, # kg/mol
    "temperature lapse rate in the troposphere": 0.0065, # K/m

    # for computation purposes
    "number of timesteps": 100, # number of timesteps in burntime (only accounts for the part of the burn where there is a thrust). Ascent timesteps typically 300-400
    "tolerated apogee difference": .01, # m, allowed space between target and real apogee
    "smallest allowed inner fuel radius": 0.01, # m; equivalent to 1 cm
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
    "po error": False, 
    "po error information": "blablabla",
}

# process inputs, convert to SI
def read_rocket_attributes(input_file):
    rocket_inputs = input_file
    # convert to sane (SI unit) values
    rocket_inputs["chamber pressure"] = float(rocket_inputs["chamber pressure"]) * 6894.76 # psi to Pa
    rocket_inputs["oxidizer mass flow rate"] = float(rocket_inputs["oxidizer mass flow rate"])
    rocket_inputs["fuel external diameter"] = float(rocket_inputs["fuel external diameter"]) * 0.0254 # convert from in to m
    rocket_inputs["fuel length"] = float(rocket_inputs["fuel length"]) * 0.0254 # convert from in to m
    rocket_inputs["target apogee"] = float(rocket_inputs["target apogee"]) * 0.3048 # ft to m
    rocket_inputs["launch site altitude"] = float(rocket_inputs["launch site altitude"]) * 0.3048 # ft to m
    rocket_inputs["dry mass"] = float(rocket_inputs["dry mass"])
    rocket_inputs["rocket external diameter"] = float(rocket_inputs["rocket external diameter"]) * 0.0254 # convert from in to m
    rocket_inputs["drag coefficient"] = float(rocket_inputs["drag coefficient"])
    rocket_inputs["launch angle"] = float(rocket_inputs["launch angle"])
    rocket_inputs["fuel grain density"] = float(rocket_inputs["fuel grain density"])
    rocket_inputs["regression rate scaling coefficient"] = float(rocket_inputs["regression rate scaling coefficient"])
    rocket_inputs["regression rate exponent"] = float(rocket_inputs["regression rate exponent"])

    return rocket_inputs
