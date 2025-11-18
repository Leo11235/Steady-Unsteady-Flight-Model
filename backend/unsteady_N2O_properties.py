'''
Full guide of variable names: 
    - p:             saturated pressure of N2O [Pa]
    - v_v:           vapor molar volume [m^3/kmol]
    - u_l:           #### liquid internal energy [kJ/(mol-K) or kJ/(mol)]
    - u_v:           #### vapor internal energy [kJ/(mol-K) or kJ/(mol)] 
    - h_l:           liquid molar enthalpy [kJ/mol]
    - h_v:           vapor molar enthalpy [kJ/mol]
    - c_v_v:         heat capacity at constant volume of vapor [kJ/(mol*K)]
    - c_p_v:         heat capacity at constant pressure of vapor [kJ/(mol*K)]
    - d_v_v/d_T
    - d_u_l/d_T
    - d_h_l/d_T
    - d_u_v/d_T
    - d_h_v/d_T
    - d_c_v_v/d_T
    - d_c_p_v/d_T
    
    - d_p/d_T # idk where this one comes from
'''


import re, json, math, bisect

# returns N2O properties dict
def initialize_N2O_properties_dict():
    # get N2O properties file
    file = "./data/N2O_perperties_lookup_table.jsonc"
    with open(file, 'r') as f:
        content = f.read()
    # remove comments
    cleaned = re.sub(r'//.*', '', content)
    cleaned = re.sub(r'/\*.*?\*/', '', cleaned, flags=re.DOTALL)
    # parse cleaned file into dicr
    N2O_properties_dict = json.loads(cleaned)
    
    # conversions (I hope these are right)
    # some quantities are in [.../(kmol)], need to be converted to [.../mol] for SI units
    N2O_properties_dict["v_v"] = [x / 1000 for x in N2O_properties_dict["v_v"]]
    N2O_properties_dict["d_v_v/d_T"] = [x / 1000 for x in N2O_properties_dict["d_v_v/d_T"]]
    N2O_properties_dict["d_u_l/d_T"] = [x / 1000 for x in N2O_properties_dict["d_u_l/d_T"]]
    N2O_properties_dict["d_h_l/d_T"] = [x / 1000 for x in N2O_properties_dict["d_h_l/d_T"]]
    N2O_properties_dict["d_u_v/d_T"] = [x / 1000 for x in N2O_properties_dict["d_u_v/d_T"]]
    N2O_properties_dict["d_h_v/d_T"] = [x / 1000 for x in N2O_properties_dict["d_h_v/d_T"]]
    N2O_properties_dict["d_c_v_v/d_T"] = [x / 1000 for x in N2O_properties_dict["d_c_v_v/d_T"]]
    N2O_properties_dict["d_c_p_v/d_T"] = [x / 1000 for x in N2O_properties_dict["d_c_p_v/d_T"]]
    
    # some quantities are in [kJ], need to be converted to [J]
    # commented some of these out bc I think they are wrong, this will be the first thing to check when things inevitably don't work perfectly
    #N2O_properties_dict["u_l"] = [x / 1000 for x in N2O_properties_dict["u_l"]]
    N2O_properties_dict["h_l"] = [x / 1000 for x in N2O_properties_dict["h_l"]]
    #N2O_properties_dict["u_v"] = [x / 1000 for x in N2O_properties_dict["u_v"]]
    N2O_properties_dict["h_v"] = [x / 1000 for x in N2O_properties_dict["h_v"]]
    N2O_properties_dict["c_v_v"] = [x / 1000 for x in N2O_properties_dict["c_v_v"]]
    N2O_properties_dict["c_p_v"] = [x / 1000 for x in N2O_properties_dict["c_p_v"]]
    
    return N2O_properties_dict


# returns property_name interpolated at tank_temp for any property_name in N2O_properties_dict
def get_N2O_property(property_name, tank_temp, N2O_properties_dict):
    # throw error if tank_temp is outside the data range
    if tank_temp < 183 or tank_temp > 309:
        raise ValueError("Tank temperature out of valid range (183 K to 309 K)")
    
    # some variables have functions rather than relying on the lookup table
    if property_name == 'p':
        return p_sat(tank_temp) # saturated pressure of N2O, Pa
    elif property_name == 'd_p/d_T':
        return dp_sat_dT(tank_temp) # saturated pressure of N2O with respect to tank temp, dPa/dT
    elif property_name == 'v_l':
        return v_l_sat(tank_temp) # saturated liquid molar volume of nitrous, m^3/mol
    elif property_name == 'd_v_l/d_T':
        return dv_l_sat_dT(tank_temp) # saturated liquid molar volume of nitrous with respect to tank temp, d(m^3/mol)/dT
    
    T_list = N2O_properties_dict['T']
    
    # for given tank_temp, find the two 'neighboring' points
    # eg if tank_temp = 196, low_index = 195 and high_index = 200
    # _position --> position of item within the list
    # _value --> value of item within that position
    low_index_position = bisect.bisect_left(T_list, tank_temp)-1
    low_index_value = T_list[low_index_position]
    high_index_position = bisect.bisect_left(T_list, tank_temp)
    high_index_value = T_list[high_index_position]
    
    # find how far tank_temp is from either end of the indices
    low_index_diff = abs(low_index_value - tank_temp)
    high_index_diff = abs(high_index_value - tank_temp)
    
    # if either index = 0, simply evaluate property_name at that index
    if high_index_diff == 0:
        return N2O_properties_dict[property_name][high_index_position]
    elif low_index_diff == 0:
        raise ValueError("low_index_diff == 0 -- shouldn't happen but if it does, just use the same logic as if (high_index_diff == 0) above")
    else:
        # interpolate best match for input tank_temp
        # match will be somewhere between low_index and high_index
        print(f'matching low index for T = {tank_temp}: {N2O_properties_dict[property_name][low_index_position]};;; diff: {low_index_diff}')
        print(f'matching high index for T = {tank_temp}: {N2O_properties_dict[property_name][high_index_position]};;; diff: {high_index_diff}')
        return (high_index_diff * N2O_properties_dict[property_name][low_index_position] + low_index_diff * N2O_properties_dict[property_name][high_index_position]) / (high_index_diff + low_index_diff)


# functions below are in a special category, they do not use the lookup table
# return saturated pressure of N2O at tank temperature T_T
def p_sat(T_T):
    # from literature
    c_1 = 96.512
    c_2 = -4045
    c_3 = -12.277
    c_4 = 0.00002886
    c_5 = 2
    
    exponent = c_1 + c_2 / T_T + c_3 * math.log(T_T) + c_4 * T_T ** c_5
    return math.exp(exponent)

# return rate of change of saturated pressure of N2O at tank temperature T_T
def dp_sat_dT(T_T):
    c_2 = -4045
    c_3 = -12.277
    c_4 = 0.00002886
    c_5 = 2
    
    return p_sat(T_T) * (-c_2/T_T**2 + c_3/T_T + c_4*c_5*T_T**(c_5-1))

# return saturated liquid molar volume of nitrous oxide at temperature T_T
def v_l_sat(T_T):
    c_1 = 2.781
    c_2 = 0.27244
    c_3 = 309.57
    c_4 = 0.2882
    
    term_1 = 1 + (1 - T_T/c_3) ** c_4
    term_2 = c_2 ** term_1
    return term_2 / (c_1 * 1000) # convert kmol to mol

# return RoC of saturated liquid molar volume of nitrous oxide at temperature T_T
def dv_l_sat_dT(T_T):
    c_2 = 0.27244
    c_3 = 309.57
    c_4 = 0.2882
    
    term_1 = (c_4 / c_3) * math.log(c_2) * v_l_sat(T_T)
    term_2 = (1 - T_T/c_3) ** (c_4 - 1)
    return term_1 * term_2 / 1000 # convert kmol to mol