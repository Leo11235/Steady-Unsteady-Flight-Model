from math import e
from pathlib import Path
import json5


_ROOT_DIR = Path(__file__).resolve().parents[3]
# example: to get to natural_constants.jsonc:
# natural_constants_file_path = _ROOT_DIR / "src" / "backend" / "static_data" / "natural_constants.jsonc"
# with open(natural_constants_file_path, "r") as f:
#     constants_dict = json5.load(f)
# print(constants_dict)

# calculates air density for a given altitude
def calculate_air_density(constants_dict, height):
    T0 = constants_dict["sea level temperature"]
    T11 = constants_dict["stratosphere temperature"]
    L = constants_dict["temperature lapse rate in the troposphere"]
    p0 = constants_dict["sea level air density"]
    g = constants_dict["sea level gravity"]
    M = constants_dict["air molar mass"]
    Ru = constants_dict["universal gas constant"]
    if height <= 11000: # if the rocket is in the troposphere
        p_air = p0 * ((T0 - L * height) / T0) ** (g*M/(Ru*L) - 1)
    else: # the rocket is in the stratosphere
        p11 = p0 * (T11 / T0) ** (g*M/(Ru*L) - 1) # pressure at 11000 m
        p_air = p11 * e ** ((11000 - height)*(g*M/(Ru*T11)))

    return p_air

# returns ambient pressure at rocket_height
def calculate_ambient_pressure(constants_dict, rocket_height):
    return

# calculates post-burn trajectory of the rocket
def calculate_flight(state_vector:dict, rocket_inputs:dict, constants_dict:dict):
    dt = rocket_inputs["timestep length"]
    sy_R = state_vector["sy_R"][-1]
    sx_R = state_vector["sx_R"][-1]
    vy_R = state_vector["vy_R"][-1]
    vx_R = state_vector["vx_R"][-1]
    ay_R = state_vector["ay_R"][-1]
    ax_R = state_vector["ax_R"][-1]

    C_d_drogue = rocket_inputs["drogue parachute drag coefficient"]
    A_drogue = rocket_inputs["drogue parachute frontal area"]
    H_deployment = rocket_inputs["main parachute deployment altitude"]
    C_d_main = rocket_inputs["main parachute drag coefficient"]
    A_main = rocket_inputs["main parachute frontal area"]

    (F_y, F_x) = 0
    m_R = 0

    # ascent
    while(vy_R >= 0):
        #new_ay_R = F_Y

        state_vector["sy_R"].append()
        state_vector["sx_R"].append()
        state_vector["vy_R"].append()
        state_vector["vx_R"].append()
        state_vector["ay_R"].append()
        state_vector["ax_R"].append()

    # drogue chute descent
    while(vy_R >= H_deployment):

        state_vector["sy_R"].append()
        state_vector["sx_R"].append()
        state_vector["vy_R"].append()
        state_vector["vx_R"].append()
        state_vector["ay_R"].append()
        state_vector["ax_R"].append()

    # main chute descent
    while(sy_R <= rocket_inputs["launch site altitude"]):

        state_vector["sy_R"].append()
        state_vector["sx_R"].append()
        state_vector["vy_R"].append()
        state_vector["vx_R"].append()
        state_vector["ay_R"].append()
        state_vector["ax_R"].append()
    return