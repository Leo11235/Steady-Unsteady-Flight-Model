import json, re, numpy, math
from unsteady_N2O_properties import get_N2O_property

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


# initialize state vector using either ullage factor or tank internal length as an input
def initialize_state_vector(rocket_inputs, N2O_properties_dict):
    # INITIALIZE CV1: tank state variables
    # initialize saturated N2O properties
    T_T_0 = rocket_inputs['tank initial temperature']
    v_l = get_N2O_property(T_T_0, 'v_l', N2O_properties_dict) # v_l = liquid molar volume at T_T_0
    v_v = get_N2O_property(T_T_0, 'v_v', N2O_properties_dict) # v_v = vapor molar volume at T_T_0
    p_0 = get_N2O_property(T_T_0, 'p', N2O_properties_dict) # pressure?
        
    # 
    #
    # METHOD OF PICKING ULLAGE OR TANK LENGTH SETUP LIKELY TO CHANGE
    #
    # unkown variables: 
        # V_l (liquid volume), 
        # n_l (tank oxidizer liquid moles)
        # n_v (tank oxidizer vapor moles)
        # L_T (tank length)
        # L_dt (dip tube length)
    #
    if rocket_inputs.get("tank internal length", "") != "":
        return initialize_state_vector_using_tank_length(rocket_inputs, T_T_0, v_l, v_v, p_0)
    elif rocket_inputs.get("tank ullage factor", "") != "":
        V_l, n_l, n_v, L_T, L_dt = initialize_state_vector_using_ullage(rocket_inputs, T_T_0, v_l, v_v, p_0)
    else:
        raise ImportError("Tank length and ullage not found")
    
    # INITIALIZE CV2: combustion chamber variables
    # unpack more variables
    R_f = rocket_inputs["fuel external radius"]
    m_f_tot = rocket_inputs["total initial fuel mass"]
    p_f = rocket_inputs["fuel density"]
    L_f = rocket_inputs["fuel length"]
    
    # calculate unknown variables
    r_f = math.sqrt(R_f**2 - m_f_tot/(math.pi*p_f*L_f))
    m_f = 0 # fuel in the chamber
    m_o = 0 # oxidizer in the chamber
    p_C = 11111111111111111111111111111 # initial chamber pressure to be calculated as a function of initial rocket height # what am I smoking it should be an input
    
    # INITIALIZE CV4: rocket (body) variables
    z_R = rocket_inputs["launch site altitude"]
    v_R = 0
    a_R = 0
    
    # first 0 is time
    return (0, n_v, n_l, T_T_0, r_f, m_o, m_f, p_C, z_R, v_R, a_R)

# uses tank length to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_tank_length(rocket_inputs, T_T_0, v_l, v_v, p_0):
    # not yet done but will use the same logic as initialize_state_vector_using_ullage()
    return

# uses tank ullage factor to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_ullage(rocket_inputs, T_T_0, v_l, v_v, p_0):
    # known variables
    m_o_tot_0 = rocket_inputs["oxidizer total initial mass"]
    W_o = rocket_inputs["oxidizer molar mass"]
    d_T = rocket_inputs["tank internal diameter"]
    D_dt = rocket_inputs["dip tube external diameter"]
    d_dt = rocket_inputs["dip tube internal diameter"]
    U = rocket_inputs["tank ullage factor"]
    
    # we need to solve for x in the A*x=b system below
    A = numpy.array([
        [0,1,1,0,0],
        [-1,v_l,0,0,0],
        [-U,0,v_v,0,0],
        [-4*U/math.pi,0,0,0,d_T**2-D_dt**2+d_dt**2],
        [-4/math.pi,0,0,d_T**2,-d_T**2]
    ])
    b = numpy.array([m_o_tot_0/W_o,0,0,0,0])
    
    V_l, n_l, n_v, L_T, L_dt = numpy.linalg.solve(A, b)
    
    print(V_l, n_l, n_v, L_T, L_dt)
    
    return V_l, n_l, n_v, L_T, L_dt