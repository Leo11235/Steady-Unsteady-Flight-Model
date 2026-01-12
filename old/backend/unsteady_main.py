import unsteady_N2O_properties, unsteady_variable_initialization, unsteady_ODE_master

# whole program is run from here
def main(input_json):
    # initialize rocket inputs file
    rocket_inputs = unsteady_variable_initialization.read_input_file(input_json) # returns rocket_inputs dict
    # initialize N2O properties dict and constants dict only once for the whole program
    N2O_properties_dict = unsteady_N2O_properties.initialize_N2O_properties_dict()
    constants_dict = unsteady_variable_initialization.initialize_natural_constants_dict()
    # initialize state vector (calculate value of all variables at t=0)
    state_vector = unsteady_variable_initialization.initialize_state_vector(rocket_inputs, N2O_properties_dict, constants_dict)

    print(f"initial state vector: {state_vector}")
        
    # run ODE master
    cached_data = unsteady_ODE_master.setup_ODE_master(rocket_inputs, constants_dict, state_vector)
    state_vector = unsteady_ODE_master.ODE_master(cached_data)

    

if __name__ == '__main__':
    # get file and run simulation
    file = "./unsteady_input_files/unsteady_input_1.jsonc"
    data = main(file)