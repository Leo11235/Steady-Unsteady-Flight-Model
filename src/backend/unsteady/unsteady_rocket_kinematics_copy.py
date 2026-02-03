from math import e
from math import cos
from math import sin
from pathlib import Path
import json5
from unsteady_variable_initialization import initialize_natural_constants_dict, read_input_file
import matplotlib.pyplot as plt


_ROOT_DIR = Path(__file__).resolve().parents[3]
# example: to get to natural_constants.jsonc:
# natural_constants_file_path = _ROOT_DIR / "src" / "backend" / "static_data" / "natural_constants.jsonc"
# with open(natural_constants_file_path, "r") as f:
#     constants_dict = json5.load(f)
# print(constants_dict)

constants_dict = initialize_natural_constants_dict()

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
def calculate_flight(state_vector:dict, rocket_inputs:dict, constants_dict:dict, force, m_R, theta):
    dt = rocket_inputs["timestep length"]
    sy_R = state_vector["sy_R"][-1]
    sx_R = state_vector["sx_R"][-1]
    vy_R = state_vector["vy_R"][-1]
    vx_R = state_vector["vx_R"][-1]
    ay_R = state_vector["ay_R"][-1]
    ax_R = state_vector["ax_R"][-1]

    C_d_rocket = rocket_inputs["rocket drag coefficient"]
    A_rocket = rocket_inputs["rocket frontal area"]
    C_d_drogue = rocket_inputs["drogue parachute drag coefficient"]
    A_drogue = rocket_inputs["drogue parachute frontal area"]
    H_deployment = rocket_inputs["main parachute deployment altitude"]
    C_d_main = rocket_inputs["main parachute drag coefficient"]
    A_main = rocket_inputs["main parachute frontal area"]
    F_g = constants_dict["sea level gravity"]

    # ascent
    while(vy_R >= 0):
        #print(vy_R)
        new_a_R = (1/m_R)*(force - (1/2)*C_d_rocket*A_rocket*calculate_air_density(constants_dict, sy_R)*(vy_R**2 + vx_R**2) - F_g*cos(theta))
        new_ay_R = new_a_R * cos(theta)
        new_ax_R = new_a_R * sin(theta)

        new_vy_R = vy_R + ay_R*dt
        new_vx_R = vx_R + ax_R*dt

        new_sy_R = sy_R + vy_R*dt
        new_sx_R = sx_R + vx_R*dt

        state_vector["sy_R"].append(new_sy_R)
        state_vector["sx_R"].append(new_sx_R)
        state_vector["vy_R"].append(new_vy_R)
        state_vector["vx_R"].append(new_vx_R)
        state_vector["ay_R"].append(new_ay_R)
        state_vector["ax_R"].append(new_ax_R)
        
        sy_R = state_vector["sy_R"][-1]
        sx_R = state_vector["sx_R"][-1]
        vy_R = state_vector["vy_R"][-1]
        vx_R = state_vector["vx_R"][-1]
        ay_R = state_vector["ay_R"][-1]
        ax_R = state_vector["ax_R"][-1]

    # drogue chute descent
    while(vy_R >= H_deployment):
        new_a_R = (1/m_R)(-F_g + (1/2)*(C_d_rocket*A_rocket + C_d_drogue*A_drogue)*calculate_air_density(constants_dict, sy_R)*(vy_R**2 + vx_R**2))
        new_ay_R = new_a_R * cos(theta)
        new_ax_R = new_a_R * sin(theta)

        new_vy_R = vy_R + ay_R*dt
        new_vx_R = vx_R + ax_R*dt

        new_sy_R = sy_R + vy_R*dt
        new_sx_R = sx_R + vx_R*dt

        state_vector["sy_R"].append(new_sy_R)
        state_vector["sx_R"].append(new_sx_R)
        state_vector["vy_R"].append(new_vy_R)
        state_vector["vx_R"].append(new_vx_R)
        state_vector["ay_R"].append(new_ay_R)
        state_vector["ax_R"].append(new_ax_R)
        
        sy_R = state_vector["sy_R"][-1]
        sx_R = state_vector["sx_R"][-1]
        vy_R = state_vector["vy_R"][-1]
        vx_R = state_vector["vx_R"][-1]
        ay_R = state_vector["ay_R"][-1]
        ax_R = state_vector["ax_R"][-1]

    # main chute descent
    while(sy_R > rocket_inputs["launch site altitude"]):
        new_a_R = (1/m_R)(-F_g + (1/2)*(C_d_rocket*A_rocket + C_d_main*A_main)*calculate_air_density(constants_dict, sy_R)*(vy_R**2 + vx_R**2))
        new_ay_R = new_a_R * cos(theta)
        new_ax_R = new_a_R * sin(theta)

        new_vy_R = vy_R + ay_R*dt
        new_vx_R = vx_R + ax_R*dt

        new_sy_R = sy_R + vy_R*dt
        new_sx_R = sx_R + vx_R*dt

        state_vector["sy_R"].append(new_sy_R)
        state_vector["sx_R"].append(new_sx_R)
        state_vector["vy_R"].append(new_vy_R)
        state_vector["vx_R"].append(new_vx_R)
        state_vector["ay_R"].append(new_ay_R)
        state_vector["ax_R"].append(new_ax_R)
        
        sy_R = state_vector["sy_R"][-1]
        sx_R = state_vector["sx_R"][-1]
        vy_R = state_vector["vy_R"][-1]
        vx_R = state_vector["vx_R"][-1]
        ay_R = state_vector["ay_R"][-1]
        ax_R = state_vector["ax_R"][-1]
    return state_vector


state_vector_test={
    "sy_R": [0],
    "sx_R": [0],
    "vy_R": [100],
    "vx_R": [10],
    "ay_R": [0],
    "ax_R": [0]
}
print("\n\n\n\n")
print(f'{_ROOT_DIR} / user_data /  simulation_configs / unsteady_input_test_1.jsonc')
rocket_inputs_test = {
    "rocket name": "Please please pleeeeease work",
    "timestep length": 0.01, # [s]
    "rocket dry mass": 90, # [kg]
    "rocket drag coefficient": 0.64, # [.]
    "rocket frontal area": 0.1257, # [m^2]
    "rocket launch angle": 6, # [°]
    "launch site altitude": 360, # z_R [m]
    "drogue parachute drag coefficient": 1.2, # [.]
    "drogue parachute frontal area": 0.28, # [m^2]
    "main parachute deployment altitude": 600, # [m]
    "main parachute drag coefficient": 1.5, # [.]
    "main parachute frontal area": 9.62, # m^2
}

state_vector_test = calculate_flight(state_vector_test, rocket_inputs_test, constants_dict, 0, 1, 0.1)
plt.plot(state_vector_test["sx_R"], state_vector_test["sy_R"])
plt.show()