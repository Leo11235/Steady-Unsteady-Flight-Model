import unsteady_N2O_properties, unsteady_variable_initialization, unsteady_ODE_master

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
    'p_C': [], # combustion chamber pressure
    # CV4: entire rocket (might split into horizontal and vertical components)
    'z_R': [], # rocket altitude
    'v_R': [], # rocket total velocity
    'a_R': [], # vertical acceleration
}

state_derivative_vector = {} # tracks the derivatives of the same variables as in the state vector

# whole program is run from here
def main(input_json):
    # initialize rocket inputs file
    rocket_inputs = unsteady_variable_initialization.read_input_file(input_json) # returns rocket_inputs dict
    # initialize N2O properties dict and constants dict only once for the whole program
    N2O_properties_dict = unsteady_N2O_properties.initialize_N2O_properties_dict()
    constants_dict = unsteady_variable_initialization.initialize_natural_constants_dict()
    # initialize state vector (calculate value of all variables at t=0)
    x_0 = unsteady_variable_initialization.initialize_state_vector(rocket_inputs, N2O_properties_dict, constants_dict)
    unsteady_ODE_master.update_state_vector(x_0)

    print(f"initial state vector: {x_0}")
    print()
    print(state_vector)
    
    # run ODE master
    cached_data = unsteady_ODE_master.setup_ODE_master(rocket_inputs, constants_dict)
    state_vector = unsteady_ODE_master.ODE_master(cached_data, state_vector)

    

if __name__ == '__main__':
    # get file and run simulation
    file = "./unsteady_input_files/unsteady_input_1.jsonc"
    data = main(file)