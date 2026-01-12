import pandas as pd 
import numpy as np
from pathlib import Path

_ROOT_DIR = Path(__file__).resolve().parent.parent.parent
# load lookup table
def load_lookup_table():
    json_path = _ROOT_DIR / "Backend" / "Static_data" / "pypropep_lookup_table.json"
    return pd.read_json(json_path, orient="table")
ppp_df_loaded = load_lookup_table()

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

def get_chamber_properties_with_partials(OF, p_C, delta_OF=0.01, delta_p=1):
    # get base values
    T_base, W_base, gamma_base = pyPROPEP_interpolation_lookup(OF, p_C, lookup_table=ppp_df_loaded)
    
    # Perturb OF (for ∂/∂(OF)) - keep pressure constant
    T_OF_plus, W_OF_plus, _ = pyPROPEP_interpolation_lookup(OF + delta_OF, p_C, lookup_table=ppp_df_loaded)
    T_OF_minus, W_OF_minus, _ = pyPROPEP_interpolation_lookup(OF - delta_OF, p_C, lookup_table=ppp_df_loaded)
    
    # Perturb pressure (for ∂/∂p_C) - keep OF constant  
    T_p_plus, W_p_plus, _ = pyPROPEP_interpolation_lookup(OF, p_C + delta_p, lookup_table=ppp_df_loaded)
    T_p_minus, W_p_minus, _ = pyPROPEP_interpolation_lookup(OF, p_C - delta_p, lookup_table=ppp_df_loaded)
    
    # Calculate partial derivatives (central difference)
    dT_dOF = (T_OF_plus - T_OF_minus) / (2 * delta_OF)
    dW_dOF = (W_OF_plus - W_OF_minus) / (2 * delta_OF)
    
    dT_dp = (T_p_plus - T_p_minus) / (2 * delta_p)
    dW_dp = (W_p_plus - W_p_minus) / (2 * delta_p)
    
    return T_base, W_base, gamma_base, dT_dOF, dW_dOF, dT_dp, dW_dp