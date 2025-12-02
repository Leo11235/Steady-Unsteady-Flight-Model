from unsteady_N2O_properties import get_N2O_property
import numpy as np
from scipy.optimize import fsolve
from PROPEP_lookup_table.pyPROPEP_simple import get_chamber_properties_with_partials, pyPROPEP_interpolation_lookup

# do all the computationally expensive one-time setup stuff here, before passing it off to ODE_master
def setup_ODE_master(rocket_inputs, constants_dict):
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
def ODE_master(cached_data, state_vector): 
    # ========================================
    # ========== Unpack cached data ==========
    # ========================================
    
    (x, dicts, dt, v_l, v_v, p_T, p_loss, C_i, N_i, A_i, W_o, a, n, p_f, L_f, A_t, V_pre, V_post) = cached_data # there will be more stuff in the paratheses
    
    # x is the state vector, also needs to be unpacked (need to verify units)
        # n_v = moles of N2O vapor in the tank [mol]
        # n_l = moles of N2O liquid in the tank [mol]
        # T_T = tank temperature [K]
        # r_f = fuel cell internal radius [m]
        # m_o = oxidizer mass in the combustion chamber [kg]
        # m_f = fuel mass in the combustion chamber [kg]
        # p_C = combustion chamber pressure [Pa]
        # z_R = rocket altitude [m]
        # v_R = rocket velocity [m/s]
        # a_R = rocket acceleration [m/s^2]
    (time, n_v, n_l, T_T, r_f, m_o, m_f, p_C, z_R, v_R, a_R) = x
    
    # other variables
        # dt: timestep length
        # v_l: Oxidizer liquid molar volume
        # v_v: Oxidizer vapor molar volume
        # p_T
        # p_loss
        # C_i, N_i, A_i: injector flow coefficient, number of holes, and area per hole
        # a: regresion rate scaling constant
        # W_o: Oxidizer molar weight
        # n: Regression rate exponent
        # p_f: Fuel density
        # L_f: fuel cell length
        # A_t: nozzle throat area
        # V_pre: pre-chamber volume
        # V_post: post-chamber volume
        
    # unpack dicts; these serve as lookup tables
    (N2O_properties_dict, PROPEP_lookup_dict) = dicts
    
    # create a time derivative state vector, with values to be filled in during the following calculations
    dx_dt = [time, 0,0,0,0,0,0,0,0,0,0]
    
    
    
    
    ########################
    # ODE while loop:
    while x[9]>0 and time != 1: # ends when velocity is 0 AND more than 1 second has elapsed
    
        # ============================================
        # ========== CV1: Tank calculations ==========
        # ============================================

        # in this section, the following quantities are calculated:
        # x = [n_v, n_l, T_T]

        # For the tank, there are two cases: liquid blowdown (n_l>0; there is liquid in the tank) and gaseous blowdown (n_l=0; there is only vapor in the tank)
        if n_l > 1e-10: # Liquid-vapor blowdown; pick a number slightly above 0 so the simulation doesn't stall
            # get N2O properties from lookup table
            d_v_v_d_T_T = get_N2O_property("d_v_v/d_T", T_T, N2O_properties_dict) # vapor molar volume with respect to temperature [m^3/(mol-K)]
            d_v_l_d_T_T = get_N2O_property("d_v_l/d_T", T_T, N2O_properties_dict) # liquid molar volume with respect to temperature [m^3/(mol-K)]
            u_v = get_N2O_property("u_v", T_T, N2O_properties_dict) # vapor internal energy [J/mol]
            u_l = get_N2O_property("u_l", T_T, N2O_properties_dict) # liquid internal energy [J/(mol-K)]
            d_u_l_d_T_T = get_N2O_property("d_u_l/d_T", T_T, N2O_properties_dict) # liquid internal energy with respect to temperature [J/(mol-K)]
            d_u_v_d_T_T = get_N2O_property("d_u_v/d_T", T_T, N2O_properties_dict) # vapor internal energy with respect to temperature [J/(mol-K)]
            h_o = get_N2O_property("h_l", T_T, N2O_properties_dict) # Oxidizer molar enthalpy (for liquid oxidizer)

            # oxidizer molar flow rate (for liquid blowdown only)
            n_dot = C_i * N_i * A_i * np.sqrt(2(p_T - p_loss - p_C) / (W_o * v_l))

            # set up Ax=b system, where dx/dt = [dn_v/dt, dn_l/dt, dT_T/dt]
            A = np.array([
                [1, 1, 0], 
                [v_v, v_l, n_v * d_v_v_d_T_T + n_l * d_v_l_d_T_T],
                [u_v, u_l, n_v * d_u_v_d_T_T + n_l * d_u_l_d_T_T]
            ])
            b = np.array([-n_dot, 0, -n_dot * h_o])

            # solve and assign values to the time derivative state vector
            x = np.linalg.solve(A, b)
            dn_v_dt, dn_l_dt, dT_T_dt = x
            dx_dt[1] = dn_v_dt
            dx_dt[2] = dn_l_dt
            dx_dt[3] = dT_T_dt

        else: # Gaseous blowdown, solve Ax=b system but different
            # polytropic exponent, approximated as the heat capacity ratio (m = heat capacity at constant pressure of vapor / heat capacity at constant volume of vapor)
            m = get_N2O_property("c_p_v", T_T, N2O_properties_dict) / get_N2O_property("c_v_v", T_T, N2O_properties_dict)

            # solve for the quantities dx/dt = [dn_v/dt, dn_l/dt, dT_T/dt]
            dn_v_dt = C_i * N_i * A_i * np.sqrt(2(p_T - p_loss - p_C) / (W_o * v_v))
            dn_l_dt = 0
            dT_T_dt = T_T * (m-1) * dn_v_dt / n_v

            # assign values to the time derivative state vector
            dx_dt[1] = dn_v_dt
            dx_dt[2] = dn_l_dt
            dx_dt[3] = dT_T_dt


        # ==========================================================
        # ========== CV2: Combustion chamber calculations ==========
        # ==========================================================

        # in this section, the following quantities are calculated:
        # x = [r_f, m_o, m_f, p_C]

        # r_f
        if r_f > 1e-10: # might be smaller than ideal
            dr_f_dt = a * ((W_o * n_dot)/(np.pi * r_f**2)) ** n
        else:
            dr_f_dt = 0

        # m_o and m_f
        OF = m_o / m_f # oxidizer to fuel ratio
        dm_n_dt = dm_o_out_dt + dm_f_out_dt # Nozzle total propellant mass flow rate

        dm_o_in_dt = W_o * n_dot # ox flow in
        dm_o_out_dt = dm_n_dt / (1+1/OF) # ox flow out
        dm_o_dt = dm_o_in_dt - dm_o_out_dt # oxidizer mass flow rate = ox flow in - ox flow out

        dm_f_in_dt = p_f * 2 * np.pi * r_f * L_f * dr_f_dt # fuel flow in
        dm_f_out_dt = dm_n_dt / (1+OF) # fuel flow out
        dm_f_dt = dm_f_in_dt - dm_f_out_dt # fuell mass flow rate = fuel flow in - fuel flow out

        # p_C
        # get PROPEP values and PROPEP partial derivatives
        T_c, W_c, gamma, dT_dOF, dW_dOF, dT_dp, dW_dp = get_chamber_properties_with_partials(OF, p_C)
        # chamber gaseous mass storage & its derivative w.r.t. time
        m_c = m_o + m_f
        dm_c_dt = dm_o_in_dt + dm_f_in_dt - dm_n_dt
        # OF ratio derivative w.r.t. time
        dOF_dt = (1/m_f)*(dm_o_dt - OF*dm_f_dt)
        # chamber volume & its derivative w.r.t. time
        V_c = V_pre + V_post + np.pi * r_f**2 * L_f
        dV_c_dt = 2* np.pi * r_f * dr_f_dt * L_f
        # finally, actually calculate p_C
        dp_C_dt = (dm_c_dt/m_c - dV_c_dt/V_c + dOF_dt(dT_dOF/T_c + dW_dOF/W_c)) / (1/p_C - dT_dp/T_c + dW_dp/W_c)

        # assign calculated values to the time derivative state vector
        dx_dt[4] = dr_f_dt
        dx_dt[5] = dm_o_dt
        dx_dt[6] = dm_f_dt
        dx_dt[7] = dp_C_dt

        # ==============================================
        # ========== CV3: Nozzle calculations ==========
        # ==============================================

        # Get M1

        term1 = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M_1**2) #Eq 3.33
        term2 = (gamma + 1) / (2 * (gamma - 1)) #eq. 3.33

        eq_3pt33 = (1 / M_1) * term1**term2
        M_1 = fsolve(eq_3pt33, 0.5)[0] # First element of the Numpy Array

        # Get p1

        term3 = (1+((gamma-1)/2)*(M_1)**2) #Every terms in the brackets of eq. 3.34
        term4 = (gamma/(gamma - 1))
        p_1 = p_C*(term3**term4)

        # Get M2x

        term5 = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M_2x**2) #Replace M2x into Eq 3.33
        term6 = (gamma + 1) / (2 * (gamma - 1)) #Replace M2x into eq. 3.33

        eq_3pt34 = (1 / M_2x) * term5**term6
        M_2x = fsolve(eq_3pt34, 0.5)[0] #First element of the Numpy Array

        # Get p2x

        term7 = (1+((gamma-1)/2)*(M_2x)**2) #Every terms in the brackets of eq. 3.34
        term8 = (gamma/(gamma - 1))
        p_2x = p_C*(term3**term4)

        #Get p_2

        top_3pt37 = 2*gamma*((M_2x)**2)- (gamma - 1)
        bottom_3pt37 = (gamma +1)
        p_2 = p_2x ( top_3pt37/bottom_3pt37)

        # Get M_2

        top_3pt38 = 2+ (gamma -1)((M_2x)**2)
        bottom_3pt38 = 2*gamma*((M_2x)**2)-(gamma-1)
        M_2 = sqrt(top_3pt38/bottom_3pt38)

        if p_amb >= p_1:
            p_e = p_amb
            if p_amb = p_1: # Case where the flow is choked
                M_t =1
                M_e = M_1
            elif p_amb > p_1: # Case where the flow is subsonic through the nozzle
                #Get M_e

                bracket_3pt35 = (1+((gamma-1)/2)*(M_e)**2) #Every terms in the brackets of eq. 3.34
                exp_3pt35 = (gamma/(gamma - 1))
                
                eq_3pt35 = bracket_3pt35**exp_3pt35 - (p_C/p_e)
                M_e = fsolve(eq_3pt35, 0.5)

                #Get M_t

                top_3pt36 = 1+((gamma-1)/2)*(M_t**2)
                bottom_3pt36 = 1+((gamma-1)/2)*(M_e**2)
                exp_3pt36 = (gamma+1)/(2(gamma-1))
                eq_3pt36 = (M_e/M_t)*(top_3pt36/bottom_3pt36)**(exp_3pt36) - (A_t/A_e)

       if (p_1 > p_amb) and (p_amb >= p_2):
        p_e = p_amb
        M_t =1 #Flow is choked

        if p_2 = p_amb:
            M_e = M_2

        elif p_2 > p_amb:
            bracket_3pt40 = (2/(gamma-1))*(1+((gamma-1)/2)*(M_e)**2)
            exp_3pt40 = (gamma+1) /(2*(gamma-1))
            eq_3pt40 = (1/M_e)*(bracket_3pt40)**(exp_3pt40) - (A_e/A_t)

            M_e = fsolve(eq_3pt40, 0.5)
        

        #Engine trust

        #Get T_e 

        rhs_3pt43 = 1 + ((gamma-1)/2)*(M_e)**2
        T_e = T_c/rhs_3pt43

        #Get v_e

        v_e = M_e* sqrt((gamma* R_u * T_e)/W_c)

        #Get m_dot_n

        exp_3pt44 = -(gamma+1)/(2(gamma-1))

        m_dot_n = A_t*p_C*M_t* sqrt((gamma* W_c)/R_u * T_e) * (1+ (((gamma-1)/2))*(M_t)**2)**exp_3pt44

        #Trust
        F = m_dot_n * v_e + (p_e - p_amb)*A_e

        








            






        


        

                
        
        









    






        # ===================================================
        # ========== CV4: Rocket body calculations ==========
        # ===================================================
        
        # in this section, the following quantities are calculated:
        # x = [z_R, v_R, a_R]


        # ==================================
        # ==== Compute new state vector ====
        # ==================================

        x[1] = dx_dt[1] * dt
        x[2] = dx_dt[2] * dt
        x[3] = dx_dt[3] * dt
        x[4] = dx_dt[4] * dt
        x[5] = dx_dt[5] * dt
        x[6] = dx_dt[6] * dt
        x[7] = dx_dt[7] * dt
        x[8] = dx_dt[8] * dt
        x[9] = dx_dt[9] * dt
        x[10] = dx_dt[10] * dt
        
        update_state_vector(x, state_vector)

    return state_vector



# updates the state vector by applying each field in input vector x to the correct place in the state_vector dict
def update_state_vector(x, state_vector):
    state_vector["time"].append(x[0])
    state_vector["n_v"].append(x[1])
    state_vector["n_l"].append(x[2])
    state_vector["T_T"].append(x[3])
    state_vector["r_f"].append(x[4])
    state_vector["m_o"].append(x[5])
    state_vector["m_f"].append(x[6])
    state_vector["p_C"].append(x[7])
    state_vector["z_R"].append(x[8])
    state_vector["v_R"].append(x[9])
    state_vector["a_R"].append(x[10])