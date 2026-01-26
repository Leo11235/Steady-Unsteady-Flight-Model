from .unsteady_N2O_properties import get_N2O_property
from ..PROPEP.interpolate_PROPEP import get_chamber_properties_with_partials
from .unsteady_rocket_kinematics import calculate_ambient_pressure
import numpy as np
from scipy.optimize import fsolve


# do all the computationally expensive one-time setup stuff here, before passing it off to ODE_master
# caching in particular ensures we are only passing the necessary inputs along to ODE_master
def setup_ODE_master(rocket_inputs, constants_dict, state_vector):
    # Rocket (body)
    m_dry = rocket_inputs["rocket dry mass"]
    C_d = rocket_inputs["rocket drag coefficient"]
    A_R = rocket_inputs["rocket frontal area"]
    theta = rocket_inputs["rocket launch angle"]
    
    # Parachutes
    C_d_drogue = rocket_inputs["drogue parachute drag coefficient"]
    A_drogue = rocket_inputs["drogue parachute frontal area"]
    H_deployment = rocket_inputs["main parachute deployment altitude"] # might rename, listed as "-" in Joel's model
    C_d_main = rocket_inputs["main parachute drag coefficient"]
    A_main = rocket_inputs["main parachute frontal area"]
    
    # Oxidizer tank
    m_o_tot = rocket_inputs["tank oxidizer mass"]
    T_T_0 = rocket_inputs["tank temperature"]
    p_T_0 = rocket_inputs["tank pressure"]
    d_T = rocket_inputs["tank internal diameter"]
    L_T = rocket_inputs["tank internal shell length"]
    V_T = rocket_inputs["tank internal volume"]
    U = rocket_inputs["tank ullage factor"]
    
    # Dip tube
    D_dt = rocket_inputs["dip tube external diameter"]
    d_dt = rocket_inputs["dip tube internal diameter"]
    L_dt = rocket_inputs["dip tube length"]
    
    # Injector and feed
    C_i = rocket_inputs["injector discharge coefficient"]
    N_i = rocket_inputs["injector number of holes"]
    A_i = rocket_inputs["injector hole area"]
    p_loss = rocket_inputs["feed pressure loss"]
    
    # Chamber
    m_f_tot = rocket_inputs["chamber fuel mass"]
    L_f = rocket_inputs["chamber fuel length"]
    p_f = rocket_inputs["chamber fuel density"]
    R_f = rocket_inputs["chamber fuel external radius"]
    a = rocket_inputs["chamber regression rate scaling constant"]
    n = rocket_inputs["chamber regression rate exponent"]
    V_pre = rocket_inputs["pre-chamber volume"]
    V_post = rocket_inputs["post-chamber volume"]
    
    # Nozzle 
    r_t = rocket_inputs["nozzle throat radius"]
    r_e = rocket_inputs["nozzle exit radius"]
    
    # Constants
    g = constants_dict["sea level gravity"]
    R_u = constants_dict["universal gas constant"]
    W_o = constants_dict["nitrous oxide molar mass"]
    # ... more to be added
    
    # Other
    dt = rocket_inputs["timestep length"]
    
    # pack up all variables
    cached_data = (state_vector, dt, m_dry, C_d, A_R, theta, C_d_drogue, A_drogue, H_deployment, C_d_main, A_main, m_o_tot, T_T_0, p_T_0, d_T, L_T, V_T, U, D_dt, d_dt, L_dt, C_i, N_i, A_i, p_loss, m_f_tot, L_f, p_f, R_f, a, n, V_pre, V_post, r_t, r_e, g, R_u, W_o)
    
    return cached_data

# runs all the ODEs in unsteady-state as optimizedly as possible
def ODE_master(cached_data, state_vector): 
    # ========================================
    # ========== Unpack cached data ==========
    # ========================================
    
    (x, dt, m_dry, C_d, A_R, theta, C_d_drogue, A_drogue, H_deployment, C_d_main, A_main, m_o_tot, T_T_0, p_T, d_T, L_T, V_T, U, D_dt, d_dt, L_dt, C_i, N_i, A_i, p_loss, m_f_tot, L_f, p_f, R_f, a, n, V_pre, V_post, r_t, r_e, g, R_u, W_o) = cached_data
    
    # calculate nozzle throat area
    
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
    
    # create a time derivative state vector, with values to be filled in during the following calculations
    dx_dt = [time, 0,0,0,0,0,0,0,0,0,0] # the rate of change of each item in the state vector
    
    ########################
    # ODE while loop:
    while x[9]>0 and time != 1: # ends when velocity is 0 AND more than 1 second has elapsed
    
        # ============================================
        # ========== CV1: Tank calculations ==========
        # ============================================

        # in this section, the following quantities are calculated:
        # x = [n_v, n_l, T_T]
        
        # get N2O property 
        v_v = get_N2O_property("v_v", T_T)
        v_l = get_N2O_property("v_l", T_T)
        
        # For the tank, there are two cases: liquid blowdown (n_l>0; there is liquid in the tank) and gaseous blowdown (n_l=0; there is only vapor in the tank)
        if n_l > 1e-10: # Liquid-vapor blowdown; pick a number slightly above 0 so the simulation doesn't stall
            # get N2O properties from lookup table
            d_v_v_d_T_T = get_N2O_property("d_v_v/d_T", T_T) # vapor molar volume with respect to temperature [m^3/(mol-K)]
            d_v_l_d_T_T = get_N2O_property("d_v_l/d_T", T_T) # liquid molar volume with respect to temperature [m^3/(mol-K)]
            u_v = get_N2O_property("u_v", T_T) # vapor internal energy [J/mol]
            u_l = get_N2O_property("u_l", T_T) # liquid internal energy [J/(mol-K)]
            d_u_l_d_T_T = get_N2O_property("d_u_l/d_T", T_T) # liquid internal energy with respect to temperature [J/(mol-K)]
            d_u_v_d_T_T = get_N2O_property("d_u_v/d_T", T_T) # vapor internal energy with respect to temperature [J/(mol-K)]
            h_o = get_N2O_property("h_l", T_T) # Oxidizer molar enthalpy (for liquid oxidizer)

            # oxidizer molar flow rate (for liquid blowdown only)
            n_dot = C_i * N_i * A_i * np.sqrt(2(p_T - p_loss - p_C) / (W_o * v_l))

            # set up Av=b system, where dx/dt = [dn_v/dt, dn_l/dt, dT_T/dt]
            A = np.array([
                [1, 1, 0], 
                [v_v, v_l, n_v * d_v_v_d_T_T + n_l * d_v_l_d_T_T],
                [u_v, u_l, n_v * d_u_v_d_T_T + n_l * d_u_l_d_T_T]
            ])
            b = np.array([-n_dot, 0, -n_dot * h_o])

            # solve and assign values to the time derivative state vector
            v = np.linalg.solve(A, b)
            dn_v_dt, dn_l_dt, dT_T_dt = v
            dx_dt[1] = dn_v_dt
            dx_dt[2] = dn_l_dt
            dx_dt[3] = dT_T_dt

        else: # Gaseous blowdown, solve Ax=b system but different
            # polytropic exponent, approximated as the heat capacity ratio (m = heat capacity at constant pressure of vapor / heat capacity at constant volume of vapor)
            m = get_N2O_property("c_p_v", T_T) / get_N2O_property("c_v_v", T_T)

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
        dm_n_dt = dm_o_out_dt + dm_f_out_dt # Nozzle total propellant mass flow rate ########### supposed to be calculated in CV3, check if this is correct

        dm_o_in_dt = W_o * n_dot
        dm_o_out_dt = dm_n_dt / (1+1/OF)
        dm_o_dt = dm_o_in_dt - dm_o_out_dt

        dm_f_in_dt = p_f * 2 * np.pi * r_f * L_f * dr_f_dt
        dm_f_out_dt = dm_n_dt / (1+OF)
        dm_f_dt = dm_f_in_dt - dm_f_out_dt

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
        
        p_amb = calculate_ambient_pressure(z_R)
        
        # define some helper functions:
        # M_1, the subsonic exit mach when the nozzle throat is choked (required to find p_1)
        def eq_3_36(M, gamma, A_ratio):
            term1 = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M**2)
            term2 = (gamma + 1) / (2 * (gamma - 1))
            return (1 / M) * term1**term2 - A_ratio
        # solve Eqs. 3.35-3.36 for fully subsonic flow (p_amb > p_1)
        def solve_subsonic_flow(p_amb, p_C, gamma, A_t, A_e):
            def eqs(vars):
                M_t, M_e = vars
                # Eq. 3.35: p_C/p_e = [1 + (gamma-1)/2 * M_e^2]^(gamma/(gamma-1))
                eq1 = p_C/p_amb - (1 + (gamma-1)/2 * M_e**2)**(gamma/(gamma-1))

                # Eq. 3.36: A_t/A_e = (M_e/M_t)*[(1+(gamma-1)/2*M_t^2)/(1+(gamma-1)/2*M_e^2)]^((gamma+1)/(2*(gamma-1)))
                term = (1 + (gamma-1)/2 * M_t**2) / (1 + (gamma-1)/2 * M_e**2)
                eq2 = A_t/A_e - (M_e/M_t) * term**((gamma+1)/(2*(gamma-1)))

                return [eq1, eq2]
            
            # Initial guess: M_t close to 1 (near choking), M_e small
            initial_guess = [0.9, 0.2]
            return fsolve(eqs, initial_guess)
        # solves equatios 3.39-3.42 in Joel's model for shock inside the nozzle (p_2 < p_amb < p_1)
        def solve_shock_in_nozzle(p_amb, p_C, gamma, A_ratio, M_2x):
            def shock_eqs(vars):
                M_x, M_e = vars  # M_x = Mach before shock, M_e = exit Mach after shock
                # Eq. 3.39: p0y/pC = [(2γM_x^2 - (γ-1))/(γ+1)]^(-1/(γ-1)) * [(γ+1)M_x^2/(2+(γ-1)M_x^2)]^(γ/(γ-1))
                term1 = ((2*gamma*M_x**2 - (gamma-1))/(gamma+1))**(-1/(gamma-1))
                term2 = ((gamma+1)*M_x**2/(2 + (gamma-1)*M_x**2))**(gamma/(gamma-1))
                p0y_over_pC = term1 * term2

                # Eq. 3.41: p0y/p_e = [1 + (gamma-1)/2 * M_e^2]^(gamma/(gamma-1))
                p0y_over_pe = (1 + (gamma-1)/2 * M_e**2)**(gamma/(gamma-1))

                # Since p_e = p_amb and p0y/p_e = (p0y/pC) * (pC/p_amb)
                eq1 = p0y_over_pe - p0y_over_pC * (p_C/p_amb)

                # Eq. 3.40: A_e/A*_y = (1/M_e)*[(2/(γ+1))*(1+(γ-1)/2*M_e^2)]^((γ+1)/(2(γ-1)))
                # Note: A*_y is the choked throat area for the flow after the shock
                A_star_y_ratio = (1/M_e) * ((2/(gamma+1))*(1 + (gamma-1)/2*M_e**2))**((gamma+1)/(2*(gamma-1)))

                # Relationship: A_e/A*_y = (A_e/A_t) / (A_t/A*_y) where A_t/A*_y = p0y/pC
                eq2 = A_star_y_ratio - A_ratio / p0y_over_pC

                return [eq1, eq2]
            # Initial guess: shock not too far from exit
            # M_x slightly less than M_2x, M_e subsonic
            initial_guess = [M_2x * 0.8, 0.3]
            return fsolve(shock_eqs, initial_guess)

        A_e = np.pi * r_e**2 # nozzle exit area
        A_t = np.pi * r_t**2 # nozzle throat area
        A_ratio = A_e / A_t # nozzle area ratio
        
        # M_1, the first critical nozzle mach number
        M_1 = fsolve(lambda M: eq_3_36(M, gamma, A_ratio), 0.5)[0] # solve for M_1 with initial guess <1
        # M_2x, the supersonic exit mach just before the shockwave (required to find p_2)
        M_2x = fsolve(lambda M: eq_3_36(M, gamma, A_ratio), 2.0)[0] # solve for M_2x with initial guess >1
        
        # p_1, first critical pressure
        p_1 = p_C / ((1 + (gamma - 1) / 2 * M_1**2) ** (gamma / (gamma - 1)))
        # p_2x, the second critical pressure right before the shockwave
        p_2x = p_C / ((1 + (gamma - 1) / 2 * M_2x**2) ** (gamma / (gamma - 1)))
        # p_2, second critical pressure
        p_2 = p_2x * ((2 * gamma * M_2x**2 - (gamma - 1)) / (gamma + 1))
        
        # M_2
        M_2 = np.sqrt((2+(gamma-1) * M_2x**2) / (2*gamma*M_2x**2-(gamma-1)))
        
        # determine the flow regime
        if p_2 >= p_1: # first, detect physical impossibility
            raise Exception("Nozzle flow regimes not regiming (p_2 >= p_1; shouldn't happen)")
        TOLERANCE = 1e-6 # might be too small, idk
        # case 1: choked throat, subsonic diverging
        if abs(p_amb - p_1) < TOLERANCE:
            p_e = p_amb
            M_t = 1
            M_e = M_1
        # case 2: shock at exit
        elif abs(p_amb - p_2) < TOLERANCE:
            p_e = p_amb
            M_t = 1
            M_e = M_2
        # case 3: fully supersonic, underexpanded nozzle
        elif p_amb < p_2:
            p_e = p_2x
            M_t = 1 ######### NEED TO CONFIRM #########
            M_e = M_2x
        # case 4: shock inside the nozzle
        elif p_2 < p_amb < p_1:
            p_e = p_amb
            M_t = 1
            M_x, M_e = solve_shock_in_nozzle(p_amb, p_C, gamma, A_ratio, M_2x) ######### NEED TO VERIFY #########
        # case 5: fully subsonic
        elif p_amb > p_1:
            p_e = p_amb
            M_t, M_e = solve_subsonic_flow(p_amb, p_C, gamma, A_t, A_e) ######### NEED TO VERIFY #########
        # if no case (should never happen, even when the physics aren't working)
        else: 
            raise Exception("Nozzle flow regime logic error")
        
        # T_e
        T_e = T_c / (1 + ((gamma-1)/2)*M_e**2)
        # v_e
        v_e = M_e * np.sqrt((gamma* R_u * T_e)/W_c)
        # m_d derivative w.r.t. time
        dm_n_dt = A_t*p_C*M_t* np.sqrt((gamma* W_c)/(R_u * T_c)) * (1+((gamma-1)/2)*(M_t)**2)**(-(gamma+1)/(2(gamma-1)))
        # F, thrust
        F = dm_n_dt * v_e + (p_e - p_amb)*A_e
        
        
        
        # M_e, the actual mach number (required to calculate exhaust exist velocity)
        # M_t, the throat mach number (used to calculate mass flow rate)
        # v_e, exhaust exist velocity
        # dm_n_dt, propellant mass flow rate


        # ===================================================
        # ========== CV4: Rocket body calculations ==========
        # ===================================================




        # ====================================================
        # ==== Compute newest entries to the state vector ====
        # ====================================================
        
        delta_x = [dt,0,0,0,0,0,0,0,0,0,0] # how much each quantity changes in the time dt
        delta_x[1] = state_vector[-1] + dx_dt[1] * dt # n_v = moles of N2O vapor in the tank [mol]

        x[1] = dx_dt[1] * dt # n_v = moles of N2O vapor in the tank [mol]
        x[2] = dx_dt[2] * dt # n_l = moles of N2O liquid in the tank [mol]
        x[3] = dx_dt[3] * dt # T_T = tank temperature [K]
        x[4] = dx_dt[4] * dt # r_f = fuel cell internal radius [m]
        x[5] = dx_dt[5] * dt # m_o = oxidizer mass in the combustion chamber [kg]
        x[6] = dx_dt[6] * dt # m_f = fuel mass in the combustion chamber [kg]
        x[7] = dx_dt[7] * dt # p_C = combustion chamber pressure [Pa]
        # not done (rocket kinematics)
        x[8] = dx_dt[8] * dt # z_R = rocket altitude [m]
        x[9] = dx_dt[9] * dt # v_R = rocket velocity [m/s]
        x[10] = dx_dt[10] * dt # a_R = rocket acceleration [m/s^2]
        
        
        
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