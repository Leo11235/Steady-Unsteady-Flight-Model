# combustion chamber control volume: computes [r_f, m_o, m_f, p_C] for each timestep


def chamber_CV2(state_vector, rocket_inputs):
    # unpack state vector
    
    # unpack rocket_inputs values
    a = rocket_inputs["regression rate scaling coefficient"]
    n = rocket_inputs["regression rate exponent"]
    W_o = rocket_inputs["oxidizer molar weight"]
    C_i = rocket_inputs["injector discharge coefficient"]
    N_i = rocket_inputs["injector number of holes"]
    A_i = rocket_inputs["injector hole area"]
    
    
    
    compute_fuel_grain_regression_rate(a, n, W_o)




def compute_fuel_grain_regression_rate(a, n, W_o):
    return 