import unsteady_N2O_properties

state_vector = {
    'time': [],
    # CV1: tank
    'n_v': [], # moles of N2O in vapor phase the tank
    'n_l': [], # moles of N2O in liquid phase the tank
    'T_T': [], # tank temperature
    # CV2: combustion chamber
    'r_f': [], # fuel grain port internal radius
    'm_o': [], # oxidizer mass in the combustion chamber
    'm_f': [], # fuel mass in the combustion chamber
    'p_c': [], # combustion chamber pressure
    # CV4: entire rocket (might split into horizontal and vertical components)
    'z_R': [], # rocket altitude
    'v_R': [], # rocket total velocity
    'a_R': [], # vertical acceleration
}

# initialize state vector using either ullage factor or tank internal length as an input
def initialize_state_vector(rocket_inputs):
    m_o_tot = rocket_inputs["oxidizer total mass in propulsion system"]
    W_o = rocket_inputs["oxidizer molar weight"]
    U = rocket_inputs["tank ullage factor"]
    d_T = rocket_inputs["tank internal diameter"]
    D_dt = rocket_inputs["dip tube external diameter"]
    d_dt = rocket_inputs["dip tube internal diameter"]
    
    v_v = rocket_inputs["moles of N2O in vapor phase the tank"]
    v_l = rocket_inputs["moles of N2O in liquid phase the tank"]
    return



if __name__ == '__main__':
    # initialize N2O properties dict only once for the whole program
    N2O_properties_dict = unsteady_N2O_properties.initialize_N2O_properties_dict()
    
    p = unsteady_N2O_properties.get_N2O_property(186, 'v_v', N2O_properties_dict)
    print(p)