'''
Documentation: 
https://github.com/Leo11235/pypropep 
'''

import pypropep as ppp
import pandas as pd 
import numpy as np
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


#if __name__ == '__main__':
    # OF = 5, CP = 600 psi
    #print(call_PROPEP(5, 600, CP_units="psi"))
    # Chamber temp = 2954
    # Molar weight = 0.02306
    # Heat ratio = 1.227

ppp_df_loaded = pd.read_json("backend\\PROPEP_lookup_table\\pypropep_lookup_table.json", orient="table")


def pyPROPEP_interpolation_lookup(OF, CP, lookup_table=ppp_df_loaded):

    OF_str = f"{OF:.1f}"

    try:
        row = lookup_table.loc[(OF_str, CP)]
        #chamber_temp, molar_weight, heat_ratio
        return (row["chamber_temp"], row["molar_weight"], row["heat_ratio"])
    except KeyError:
        pass  #interpolate

    unique_OF = np.array(sorted({float(o) for o in lookup_table.index.get_level_values("OF")}))
    unique_CP = np.array(sorted({int(c)    for c in lookup_table.index.get_level_values("CP")}))
    
    if OF < unique_OF.min() or OF > unique_OF.max():
        raise ValueError(f"OF={OF} is outside the table bounds [{unique_OF.min()}, {unique_OF.max()}]")
    if CP < unique_CP.min() or CP > unique_CP.max():
        raise ValueError(f"CP={CP} is outside the table bounds [{unique_CP.min()}, {unique_CP.max()}]")

    OF_low  = unique_OF[unique_OF <= OF].max()
    OF_high = unique_OF[unique_OF >= OF].min()
    CP_low  = unique_CP[unique_CP <= CP].max()
    CP_high = unique_CP[unique_CP >= CP].min()

    OF_low_str  = f"{OF_low:.1f}"
    OF_high_str = f"{OF_high:.1f}"

    def get_values(OFs, CPs):
        row = lookup_table.loc[(OFs, CPs)]
        return np.array([
            row["chamber_temp"],
            row["molar_weight"],
            row["heat_ratio"]
        ], dtype=float)

    Q11 = get_values(OF_low_str,  CP_low)
    Q12 = get_values(OF_low_str,  CP_high)
    Q21 = get_values(OF_high_str, CP_low)
    Q22 = get_values(OF_high_str, CP_high)

    def bilinear_interp(x, y, x1, x2, y1, y2, Q11, Q12, Q21, Q22):
        return (Q11 * (x2 - x) * (y2 - y) + Q21 * (x - x1) * (y2 - y) + Q12 * (x2 - x) * (y - y1) + Q22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))

    result = bilinear_interp(OF, CP, OF_low, OF_high, CP_low, CP_high, Q11, Q12, Q21, Q22)

    #chamber temp, molar weight, heat ratio
    return (float(result[0]), float(result[1]), float(result[2]))

#just checking
#print(pyPROPEP_interpolation_lookup(2, 300))
#print(pyPROPEP_interpolation_lookup(3.34, 756))