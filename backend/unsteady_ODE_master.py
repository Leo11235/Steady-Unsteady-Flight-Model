from unsteady_N2O_properties import get_N2O_property
import numpy as np

# do all the computationally expensive one-time setup stuff here, before passing it off to ODE_master
def setup_ODE_master(state_vector, state_vector_history, rocket_inputs, constants_dict, N2O_properties_dict):
    # unpack rocket inputs
    # Tank
    # Chamber
    a = rocket_inputs["regression rate scaling coefficient"]
    n = rocket_inputs["regression rate exponent"]
    W_o = rocket_inputs["oxidizer molar weight"]
    C_i = rocket_inputs["injector discharge coefficient"]
    N_i = rocket_inputs["injector number of holes"]
    A_i = rocket_inputs["injector hole area"]
    # Nozzle
    # Rocket
    
    # unpack physical constants
    g = constants_dict["sea level gravity"]
    R_u = constants_dict[""]
    
    cached_data = () # all variables will be lacked into here to then be quickly unpacked
    
    return cached_data

# runs all the ODEs in unsteady-state as optimizedly as possible
def ODE_master(cached_data): 
    # ========================================
    # ========== Unpack cached data ==========
    # ========================================
    
    (x, dicts, v_l, v_v, p_T, p_loss, C_i, N_i, A_i, W_o) = cached_data # there will be more stuff in the paratheses
    
    # x is the state vector, also needs to be unpacked (need to verify units)
        # n_v = moles of N2O vapor in the tank [mol]
        # n_l = moles of N2O liquid in the tank [mol]
        # T_T = tank temperature [K]
        # r_f = fuel cell internal radius [m]
        # m_o = oxidizer mass in the combustion chamber [kg]
        # m_f = fuel mass in the combustion chamber [kg]
        # p_C = combustion chamber pressure [Pa]
        # z_R = rocket altitude
        # v_R = rocket velocity
        # a_R = rocket acceleration
    (time, n_v, n_l, T_T, r_f, m_o, m_f, p_C, z_R, v_R, a_R) = x
    
    # other variables
        # v_l = Oxidizer liquid molar volume
        # v_v = Oxidizer vapor molar volume
        
    # unpack dicts; these serve as lookup tables
    (N2O_properties_dict, PROPEP_lookup_dict) = dicts
    
    # ============================================
    # ========== CV1: Tank calculations ==========
    # ============================================
    
    # calculate some important values first
    n_dot = C_i * N_i * A_i * np.sqrt(2(p_T - p_loss - p_C) / (W_o * v_l))
    # oxidizer molar flow rate
    
    # For the tank, there are two cases: liquid blowdown (n_l>0; there is liquid in the tank) and gaseous blowdown (n_l=0; there is only vapor in the tank)
    if n_l > 1e-10: # Liquid-vapor blowdown; pick a number slightly above 0 so the simulation doesn't stall
        # get N2O properties from lookup table
        d_v_v_d_T_T = get_N2O_property("d_v_v/d_T", T_T, N2O_properties_dict) # vapor molar volume with respect to temperature [m^3/(mol-K)]
        d_v_l_d_T_T = get_N2O_property("d_v_l/d_T", T_T, N2O_properties_dict) # liquid molar volume with respect to temperature [m^3/(mol-K)]
        u_v = get_N2O_property("u_v", T_T, N2O_properties_dict) # vapor internal energy [J/mol]
        u_l = get_N2O_property("u_l", T_T, N2O_properties_dict) # liquid internal energy [J/(mol-K)]
        d_u_l_d_T_T = get_N2O_property("d_u_l/d_T", T_T, N2O_properties_dict) # liquid internal energy with respect to temperature [J/(mol-K)]
        d_u_v_d_T_T = get_N2O_property("d_u_v/d_T", T_T, N2O_properties_dict) # vapor internal energy with respect to temperature [J/(mol-K)]
        h_o = get_N2O_property # Oxidizer molar enthalpy
        
        # solve Ax=b system, where x = [n_v, n_l, T_T]
        A = np.array([
            [1, 1, 0], 
            [v_v, v_l, n_v * d_v_v_d_T_T + n_l * d_v_l_d_T_T],
            [u_v, u_l, n_v * d_u_v_d_T_T + n_l * d_u_l_d_T_T]
        ])
        b = np.array([-n_dot, 0, -n_dot * h_o])
        pass
    else: # Gaseous blowdown, solve Ax=b system but different
        pass
    
    
    
    # ==========================================================
    # ========== CV2: Combustion chamber calculations ==========
    # ==========================================================
    
    
    # ==============================================
    # ========== CV3: Nozzle calculations ==========
    # ==============================================
    
    
    # ===================================================
    # ========== CV4: Rocket body calculations ==========
    # ===================================================
    return