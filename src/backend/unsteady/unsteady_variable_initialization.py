import json, re, numpy, math
from .unsteady_N2O_properties import get_N2O_property
from .unsteady_rocket_kinematics import calculate_air_density
from pathlib import Path

_ROOT_DIR = Path(__file__).resolve().parents[3]

# process inputs, modifies simulation_settings_dict and constants_dict, and outputs a rocket_inputs dict
def read_input_file(input_file):
    # extract json from input file
    with open(input_file, 'r') as f:
        content = f.read()
    # remove comments
    content = re.sub(r'//.*?$|/\*.*?\*/', '', content, flags=re.MULTILINE | re.DOTALL)
    # parse the cleaned JSON string
    input_file =json.loads(content)
        
    # unpack input_file dicts
    #simulation_settings_inputs = input_file['simulation settings']
    #simulation_constants_inputs = input_file['simulation constants']
    rocket_inputs = input_file['rocket attributes']
        
    # override simulation settings and constants as needed        
    #simulation_settings_dict.update({k: v for k, v in simulation_settings_inputs.items() if k in simulation_settings_dict})
    #constants_dict.update({k: v for k, v in simulation_constants_inputs.items() if k in constants_dict})
    
    return rocket_inputs

# returns a dict of natural constants used throughout the simulation. 
def initialize_natural_constants_dict():
    print(_ROOT_DIR)
    # find and open file
    file = _ROOT_DIR / "src" / "backend" / "static_data" / "natural_constants.jsonc"
    with open(file, 'r') as f:
        content = f.read()
    # remove comments
    cleaned = re.sub(r'//.*', '', content)
    cleaned = re.sub(r'/\*.*?\*/', '', cleaned, flags=re.DOTALL)
    # parse cleaned file into dicr
    constants_dict = json.loads(cleaned)
    return constants_dict


# initialize state vector using either ullage factor or tank internal length as an input
def initialize_state_vector(rocket_inputs, N2O_properties_dict, constants_dict):
    # INITIALIZE CV1: tank state variables [n_v, n_l, T_T]
    # initialize saturated N2O properties
    T_T_0 = rocket_inputs['tank initial temperature']
    v_l = get_N2O_property('v_l', T_T_0, N2O_properties_dict) # v_l = liquid molar volume at T_T_0
    v_v = get_N2O_property('v_v', T_T_0, N2O_properties_dict) # v_v = vapor molar volume at T_T_0
    p_0 = get_N2O_property('p', T_T_0, N2O_properties_dict) # pressure?
    # unpack some rocket parameters
    m_o_tot_0 = rocket_inputs["oxidizer total initial mass"]
    W_o = rocket_inputs["oxidizer molar mass"]
    d_T = rocket_inputs["tank internal diameter"]
    D_dt = rocket_inputs["dip tube external diameter"]
    d_dt = rocket_inputs["dip tube internal diameter"]
    
    # decide whether to initialize tank variables using ullage or tank length
    if rocket_inputs.get("tank internal length", "") != "":
        print("Initializing state vector using tank internal length")
        return initialize_state_vector_using_tank_length(rocket_inputs, v_l, v_v, m_o_tot_0, W_o, d_T, D_dt, d_dt)
    elif rocket_inputs.get("tank ullage factor", "") != "":
        # print("Initializing state vector using tank ullage")
        V_l, n_l, n_v, L_T, L_dt = initialize_state_vector_using_ullage(rocket_inputs, v_l, v_v, m_o_tot_0, W_o, d_T, D_dt, d_dt)
    else:
        raise ImportError("Tank length and ullage not found")
    
    # INITIALIZE CV2: combustion chamber variables [r_f, m_o, m_f, p_C]
    # unpack more variables
    R_f = rocket_inputs["fuel external radius"]
    m_f_tot = rocket_inputs["total initial fuel mass"]
    p_f = rocket_inputs["fuel density"]
    L_f = rocket_inputs["fuel length"]
    
    # calculate unknown variables
    r_f = math.sqrt(R_f**2 - m_f_tot/(math.pi*p_f*L_f)) # initial fuel port internal radius
    m_f = 0 # initial fuel in the chamber
    m_o = 0 # initial oxidizer in the chamber
    p_C = calculate_air_density(constants_dict, rocket_inputs["launch site altitude"]) # initial chamber pressure; initially the same as ambient pressure
    
    # INITIALIZE CV4: rocket (body) variables [z_R, v_R, a_R]
    z_R = rocket_inputs["launch site altitude"] # height
    v_R = 0 # velocity
    a_R = 0 # acceleration
    
    return {
        'time': [0],
        # CV1: tank
        'n_v': [n_v], # moles of N2O in vapor phase the tank
        'n_l': [n_l], # moles of N2O in liquid phase the tank
        'T_T': [T_T_0], # tank temperature
        # CV2: combustion chamber
        'r_f': [r_f], # fuel cell internal radius
        'm_o': [m_o], # oxidizer mass in the combustion chamber
        'm_f': [m_f], # fuel mass in the combustion chamber
        'p_C': [p_C], # combustion chamber pressure
        # CV4: entire rocket (might split into horizontal and vertical components)
        'z_R': [z_R], # rocket altitude
        'v_R': [v_R], # rocket total velocity
        'a_R': [a_R], # vertical acceleration
    }

# uses tank ullage factor to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_ullage(rocket_inputs, v_l, v_v, m_o_tot_0, W_o, d_T, D_dt, d_dt):
    # get rocket ullage    
    U = rocket_inputs["tank ullage factor"]
    
    # need to solve for x in the A*x=b system below
    A = numpy.array([
        [0,             1,      1,       0,          0                         ],
        [-1,            v_l,    0,       0,          0                         ],
        [-U,            0,      v_v,     0,          0                         ],
        [-4*U/math.pi,  0,      0,       0,          d_T**2 - D_dt**2 + d_dt**2],
        [-4/math.pi,    0,      0,       d_T**2,    -d_T**2                    ]
    ])
    b = numpy.array([m_o_tot_0/W_o,   0,0,0,0])    
    
    # solve the system Ax=b for x
    V_l, n_l, n_v, L_T, L_dt = numpy.linalg.solve(A, b)
    
    # print()
    # print(f'ullage init vector: {V_l, n_l, n_v, L_T, L_dt}')
    # print()
    
    return V_l, n_l, n_v, L_T, L_dt

# uses tank length to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_tank_length(rocket_inputs, v_l, v_v, m_o_tot_0, W_o, d_T, D_dt, d_dt):
    # unpack rocket length
    L_T = rocket_inputs["tank internal length"]
    
    # do the lin alg stuff
    A = numpy.array([
        [0,         1,     1,      0,              0                         ], 
        [-1,        v_l,   0,      0,              0                         ], 
        [0,         0,     v_v,   -1,              0                         ], 
        [0,         0,     0,      -4/math.pi,     d_T**2 - D_dt**2 + d_dt**2], 
        [4/math.pi, 0,     0,      0,              d_T**2                    ]
    ])
    b = numpy.array([m_o_tot_0/W_o,  0,0,0, d_T**2 * L_T])
    
    V_l, n_l, n_v, V_V, L_dt = numpy.linalg.solve(A, b)
    
    return V_l, n_l, n_v, V_V, L_dt