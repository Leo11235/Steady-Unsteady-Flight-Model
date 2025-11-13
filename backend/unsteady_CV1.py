import numpy as np
from unsteady_N2O_properties import get_N2O_property
# Tank control volume: computes [n_v, n_l, T_T] for each timestep

def tank_CV1(liquid_oxidizer_volume):
    if liquid_oxidizer_volume > 0: # liquid blowdown
        return liquid_oxidizer_blowdown()
    elif liquid_oxidizer_volume == 0: # gaseous blowdown
        return gaseous_oxidizer_blowdown()
    else:
        print("Error in CV1: liquid_oxidizer_volume neither >0 or =0")

def liquid_oxidizer_blowdown(C_i, N_i, A_i,n_l, n_v, p_loss, p_c, W_o,TT, N2O_properties_dict):

    # -----------------------
    # Get properties at TT
    # -----------------------
    p_t = get_N2O_property('p', TT, N2O_properties_dict)
    v_l = get_N2O_property('v_l',TT, N2O_properties_dict)     
    v_v = get_N2O_property('v_v',TT, N2O_properties_dict)
    u_l = get_N2O_property('u_l', TT, N2O_properties_dict)       
    u_v = get_N2O_property('u_v', TT, N2O_properties_dict)
    h_o = get_N2O_property('h_o', TT, N2O_properties_dict)

    # Derivatives wrt T (these can come from spline derivatives or finite diff)
    dv_l_dT = get_N2O_property('dv_l_dT', TT, N2O_properties_dict)
    dv_v_dT = get_N2O_property('dv_v_dT', TT, N2O_properties_dict)
    du_l_dT = get_N2O_property('du_l_dT', TT, N2O_properties_dict)
    du_v_dT = get_N2O_property('du_v_dT', TT, N2O_properties_dict)

    # -----------------------
    # Eq (3.3): molar discharge rate
    # -----------------------
    
    delta_p = max(p_t - p_loss - p_c, 0.0)
    n_dot = C_i * N_i * A_i * np.sqrt(2.0 * delta_p / (W_o * v_l))

    # -----------------------
    # Build linear system A x = b for:
    # x = [dnv_dt, dnl_dt, dTT_dt]
    # -----------------------
    # (3.3) -> 1*dnv + 1*dnl + 0*dT = n_dot
    a11, a12, a13 = 1.0, 1.0, 0.0
    b1 = n_dot

    # (3.6) -> v_v*dnv + v_l*dnl + (nv*dv_v/dT + nl*dv_l/dT)*dT = 0
    a21 = v_v
    a22 = v_l
    a23 = n_v * dv_v_dT + n_l * dv_l_dT
    b2 = 0.0

    # (3.8) -> u_v*dnv + u_l*dnl + (nl*du_l/dT + nv*du_v/dT)*dT = n_dot * h_o
    a31 = u_v
    a32 = u_l
    a33 = n_l * du_l_dT + n_v * du_v_dT
    b3 = n_dot * h_o

    A = np.array([[a11, a12, a13],
                  [a21, a22, a23],
                  [a31, a32, a33]], dtype=float)

    b = np.array([b1, b2, b3], dtype=float)

    # Solve linear system
    try:
        x = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        x, *_ = np.linalg.lstsq(A, b, rcond=None)

    dnv_dt, dnl_dt, dTT_dt = x
    return dnv_dt, dnl_dt, dTT_dt

def gaseous_oxidizer_blowdown():
    return 