from .unsteady_N2O_properties import get_N2O_property
from ..PROPEP.interpolate_PROPEP import get_chamber_properties_with_partials
from .unsteady_rocket_kinematics import calculate_air_pressure
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
def ODE_master(cached_data, state_vector, constants_dict): 
    # short explanation of some otherwise confusing variables: 
        # state_vector: a dict of lists containing the various quantities tracked by the program, this is the same state vector as everywhere else in the program
        # dx_dt[i]: a list tracking the rate of change of i, where the i'th element is the same as the i'th element in state_vector.
        # dt: float, the timestep increment magnitude
    
    # ========================================
    # ========== Unpack cached data ==========
    # ========================================
    
    (state_vector, dt, m_dry, C_d, A_R, theta, C_d_drogue, A_drogue, H_deployment, C_d_main, A_main, m_o_tot, T_T_0, p_T, d_T, L_T, V_T, U, D_dt, d_dt, L_dt, C_i, N_i, A_i, p_loss, m_f_tot, L_f, p_f, R_f, a, n, V_pre, V_post, r_t, r_e, g, R_u, W_o) = cached_data
    
    # calculate nozzle throat area
    
    # x is the state vector, also needs to be unpacked (need to verify units)
        # n_v = moles of N2O vapor in the tank [mol]
        # n_l = moles of N2O liquid in the tank [mol]
        # T_T = tank temperature [K]
        # r_f = fuel cell internal radius [m]
        # m_o = oxidizer mass in the combustion chamber [kg]
        # m_f = fuel mass in the combustion chamber [kg]
        # p_C = combustion chamber pressure [Pa]
        # OF = oxidizer to fuel ratio
        # T_C = combustion chamber temperature
        # F_x = horizontal thrust component
        # F_y = vertical thrust component
        # sy_R = rocket horizontal position ASL
        # sx_R = rocket horizontal position (launchsite=0)
        # vy_R = rocket vertial velocity
        # vx_R = rocket horizontal velocity
        # ay_R = rocket vertical acceleration
        # ax_R = rocket horizontal acceleration
        # D_x = horizontal drag force
        # D_y = vertical drag force

    (time, n_v, n_l, T_T, r_f, m_o, m_f, p_C, OF, T_C, F_x, F_y, sy_R, sx_R, vy_R, vx_R, ay_R, ax_R, D_x, D_y) = state_vector
    
    # create a time derivative state vector, with values to be filled in during the following calculations
    dx_dt = [time, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # the rate of change of each item in the state vector
    
    # set some initial values
    dm_n_dt = 0 # Nozzle total propellant mass flow rate
    
    ########################
    # ODE while loop:
    #while state_vector["vy_R"][-1]<=0: # ends when velocity is less than 0
    while state_vector["n_l"][-1]>1: # temp condition for testing until the end of gaseous blowdown
        
        # First: get latest state vector values (only for the necessary items)
        n_v = state_vector["n_v"][-1]
        n_l = state_vector["n_l"][-1]
        T_T = state_vector["T_T"][-1]
        r_f = state_vector["r_f"][-1]
        m_o = state_vector["m_o"][-1]
        m_f = state_vector["m_f"][-1]
        p_C = state_vector["p_C"][-1]
        # kinematics 
        sy_R = state_vector["n_v"][-1]
        sx_R = state_vector["n_v"][-1]
        vy_R = state_vector["n_v"][-1]
        vx_R = state_vector["n_v"][-1]
        ay_R = state_vector["n_v"][-1]
        ax_R = state_vector["n_v"][-1]
        D_x = state_vector["n_v"][-1]
        D_y = state_vector["n_v"][-1]
        
        
        # disgusting assumption: if the tank temperature is far below (happens during gaseous blowdown) or above (hasn't happened in simulations), clip it
        if 185.0>T_T or T_T>305.0:
            T_T_clamped = np.clip(T_T, 185.0, 300.0)
            #print(f"clamping T_T from {T_T} to {T_T_clamped}")
            T_T = T_T_clamped
    
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
            delta_p = p_T - p_loss - p_C # pressure loss differential
            if delta_p > 0:
                n_dot = C_i * N_i * A_i * np.sqrt(2*(delta_p) / (W_o * v_l))
            else: 
                n_dot = 0.0 # prevents backflow

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
            dn_v_dt = C_i * N_i * A_i * np.sqrt(2*(p_T - p_loss - p_C) / (W_o * v_v))
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
        OF = 0 if m_f == 0 else m_o / m_f # oxidizer to fuel ratio
        #dm_n_dt = dm_o_out_dt + dm_f_out_dt # Nozzle total propellant mass flow rate ########### supposed to be calculated in CV3, check if this is correct

        dm_o_in_dt = W_o * n_dot
        dm_o_out_dt = 0 if dm_n_dt == 0 else dm_n_dt / (1+1/OF)
        dm_o_dt = dm_o_in_dt - dm_o_out_dt

        dm_f_in_dt = p_f * 2 * np.pi * r_f * L_f * dr_f_dt
        dm_f_out_dt = dm_n_dt / (1+OF)
        dm_f_dt = dm_f_in_dt - dm_f_out_dt

        # p_C
        # get PROPEP values and PROPEP partial derivatives
        T_c, W_c, gamma, dT_dOF, dW_dOF, dT_dp, dW_dp = get_chamber_properties_with_partials(OF, p_C)
        dx_dt[9] = T_c
        # chamber gaseous mass storage & its derivative w.r.t. time
        m_c = m_o + m_f
        dm_c_dt = dm_o_in_dt + dm_f_in_dt - dm_n_dt
        
        # OF ratio derivative w.r.t. time
        if m_f > 1e-8: # use a small threshold
            dOF_dt = (1/m_f) * (dm_o_dt - OF * dm_f_dt)
        else:
            dOF_dt = 0  # assume OF is constant (or forced to 1.0) while chamber is empty
        dx_dt[8] = dOF_dt
        
        # chamber volume & its derivative w.r.t. time
        V_c = V_pre + V_post + np.pi * r_f**2 * L_f
        dV_c_dt = 2* np.pi * r_f * dr_f_dt * L_f
        
        # finally, actually calculate p_C
        if m_c > 1e-7:
            dp_C_dt = (dm_c_dt/m_c - dV_c_dt/V_c + dOF_dt*(dT_dOF/T_c + dW_dOF/W_c)) / (1/p_C - dT_dp/T_c + dW_dp/W_c)
        else:
            dp_C_dt = (dm_o_in_dt + dm_f_in_dt) * (R_u / W_c) * T_c / V_c

        # assign calculated values to the time derivative state vector
        dx_dt[4] = dr_f_dt
        dx_dt[5] = dm_o_dt
        dx_dt[6] = dm_f_dt
        dx_dt[7] = dp_C_dt

        # ==============================================
        # ========== CV3: Nozzle calculations ==========
        # ==============================================
        
        p_amb = calculate_air_pressure(constants_dict, sy_R)
        
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
        dm_n_dt = A_t*p_C*M_t* np.sqrt((gamma* W_c)/(R_u * T_c)) * (1+((gamma-1)/2)*(M_t)**2)**(-(gamma+1)/(2*(gamma-1)))
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
        
        print(round(state_vector["time"][-1], 4), round(state_vector["n_l"][-1], 4))
        update_state_vector(dt, dx_dt, state_vector)

    return state_vector

# updates the state vector by applying each field in input vector x to the correct place in the state_vector dict
def update_state_vector(dt, dx_dt, state_vector): 
    # Euler timestep formula: new_x = old_x + dx_dt * dt
    state_vector["time"].append(state_vector["time"][-1] + dt) 
    # CV1: tank
    state_vector["n_v"].append(state_vector["n_v"][-1] + dx_dt[1] * dt) # n_v = moles of N2O vapor in the tank [mol]
    state_vector["n_l"].append(state_vector["n_l"][-1] + dx_dt[2] * dt) # n_l = moles of N2O liquid in the tank [mol]
    state_vector["T_T"].append(state_vector["T_T"][-1] + dx_dt[3] * dt) # T_T = tank temperature [K]
    # CV2: combustion chamber
    state_vector["r_f"].append(state_vector["r_f"][-1] + dx_dt[4] * dt) # r_f = fuel cell internal radius [m]
    state_vector["m_o"].append(state_vector["m_o"][-1] + dx_dt[5] * dt) # m_o = oxidizer mass in the combustion chamber [kg]
    state_vector["m_f"].append(state_vector["m_f"][-1] + dx_dt[6] * dt) # m_f = fuel mass in the combustion chamber [kg]
    state_vector["p_C"].append(state_vector["p_C"][-1] + dx_dt[7] * dt) # p_C = combustion chamber pressure [Pa]
    state_vector["OF"].append(state_vector["OF"][-1] + dx_dt[8] * dt) # OF = oxidizer to fuel ratio [x]
    state_vector["T_C"].append(dx_dt[9]) # T_C = combustion chamber temperature [T] /// in this one we calculate & get T_C rather than dT_C_dt
    # CV3: Nozzle
    state_vector["F_x"].append(state_vector["F_x"][-1] + dx_dt[10] * dt)
    state_vector["F_y"].append(state_vector["F_y"][-1] + dx_dt[11] * dt)
    # CV4: Rocket body
    state_vector["sy_R"].append(state_vector["sy_R"][-1] + dx_dt[12] * dt)
    state_vector["sx_R"].append(state_vector["sx_R"][-1] + dx_dt[13] * dt)
    state_vector["vy_R"].append(state_vector["vy_R"][-1] + dx_dt[14] * dt)
    state_vector["vx_R"].append(state_vector["vx_R"][-1] + dx_dt[15] * dt)
    state_vector["ay_R"].append(state_vector["ay_R"][-1] + dx_dt[16] * dt)
    state_vector["ax_R"].append(state_vector["ax_R"][-1] + dx_dt[17] * dt)
    state_vector["D_x"].append(state_vector["D_x"][-1] + dx_dt[18] * dt)
    state_vector["D_y"].append(state_vector["D_y"][-1] + dx_dt[19] * dt)


