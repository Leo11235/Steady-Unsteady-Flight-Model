'''
Documentation: 
https://github.com/Leo11235/pypropep 
'''

import pypropep as ppp
ppp.init() # initialize pypropep once for many calls

def call_PROPEP(OF, CP, CP_units = "atm"): # oxidizer to fuel ratio and chamber pressure; units can also be Pa, psi, or atm and the conversion will be handled automatically
    fuel_type = "EICOSANE (PARAFFIN)"
    oxidizer_type = "NITROUS OXIDE"
    
    if CP_units.lower() == "psi":
        CP = CP/14.696 # convert from psi to atm
    elif CP_units.lower() == "Pa":
        CP = CP*101325 # convert from Pa to atm
    elif CP_units.lower() == "atm":
        None
    else:
        raise ValueError("Provide a valid argument for CP_unit")
    
    # make sure the propellants match nicely
    try:
        fuel = ppp.PROPELLANTS[fuel_type.upper()]
        oxidizer = ppp.PROPELLANTS[oxidizer_type.upper()]
    except KeyError as e:
        return {"error": f"Invalid propellant name: {e}"}
    
    # create equilibrium object
    eq = ppp.Equilibrium()
    try: 
        # add propellants
        eq.add_propellants_by_mass([(fuel, 1.0), (oxidizer, OF)])
        # set chamber state
        eq.set_state(P=CP)
    except Exception as e:
        return {"error": f"Error setting up equilibrium state: {e}"}
    
    # get chamber temp, molar weight, and heat ratio
    try: 
        properties = eq.properties
        chamber_temp = getattr(properties, "T", None)
        molar_weight = getattr(properties, "M", None) / 1000 # convert to kg/mol
        # heat ratio is a bit more complicated
        Cp = getattr(properties, "Cp", None)
        Cv = getattr(properties, "Cv", None)
        heat_ratio = Cp / Cv if Cp is not None and Cv is not None else None
        
        return (chamber_temp, molar_weight, heat_ratio)

    except Exception as e:
        return {"error": f"pyPROPEP - Error extracting properties: {e}"}


if __name__ == '__main__':
    # OF = 5, CP = 600 psi
    print(call_PROPEP(5, 600, CP_units="psi"))
    # Chamber temp = 2954
    # Molar weight = 0.02306
    # Heat ratio = 1.227