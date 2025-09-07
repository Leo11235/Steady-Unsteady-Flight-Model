import steady_backend.fuel_mass_convergence, steady_backend.parametric_studies, steady_backend.initialize_variables, steady_backend.unit_conversion

def run_steady(rocket_inputs, parametric_inputs, simulation_type, simulation_settings):
    print(f'{simulation_settings["debug comments"]} turned on') if simulation_settings["debug comments"] != "none" else None
    # convert inputs & create a new dictionary with just the SI values for each entry (parametric inputs are converted in parametric_studies.py)
    rocket_input_values = steady_backend.initialize_variables.initialize_rocket_input_values_dict(rocket_inputs)
        
    # initialize constants dict
    constants_dict = steady_backend.initialize_variables.initialize_constants_dict(rocket_input_values["planet"])
    
    # run simulation
    if simulation_type.get("optimize values for unsteady"): # get optimal Mo, Lf, Re for unsteady use
        None
    elif simulation_type.get("simple fuel mass convergence"): # converge on fuel mass
        return steady_backend.fuel_mass_convergence.iterate_over_inner_radius(rocket_input_values, constants_dict, simulation_settings)
        # output inludes {rocket inputs, constants dict, rocket parameters, flight data}
    elif simulation_type.get("parametric study"): # parametric study with x amount of inputs
        return steady_backend.parametric_studies.run_parametric_study(rocket_input_values, constants_dict, parametric_inputs, simulation_settings)
    
    