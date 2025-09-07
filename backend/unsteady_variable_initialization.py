import json, re
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
    # initialize saturated N2O properties
    T_T_0 = float(rocket_inputs['tank initial temperature'])
    v_l = get_N2O_property(T_T_0, 'v_l', N2O_properties_dict) # v_l = liquid molar volume at T_T_0
    v_v = get_N2O_property(T_T_0, 'v_v', N2O_properties_dict) # v_v = vapor molar volume at T_T_0
    p_0 = get_N2O_property(T_T_0, 'p', N2O_properties_dict) # pressure?
    
    print(T_T_0, v_l, v_v, p_0)
    
    if rocket_inputs.get("tank internal length", "") != "":
        return initialize_state_vector_using_tank_length(rocket_inputs, T_T_0, v_l, v_v, p_0)
    elif rocket_inputs.get("tank ullage factor", "") != "":
        return initialize_state_vector_using_tank_length(rocket_inputs, T_T_0, v_l, v_v, p_0)
    else:
        raise ImportError("Tank length and ullage not found")


# uses tank length to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_tank_length(rocket_inputs, T_T_0, v_l, v_v, p_0):
    return

# uses tank ullage factor to initialize the state vector, returns state vector at t=0
def initialize_state_vector_using_ullage(rocket_inputs, T_T_0, v_l, v_v, p_0):
    
    
    return