import pandas as pd
import numpy as np
import pypropep as ppp
from pyPROPEP_simple import call_PROPEP

ppp.init() 

pypropep_data = []

for OF in np.arange(1, 10.1, 0.1):
    for CP in range(100, 1210, 10):
        chamber_temp, molar_weight, heat_ratio = call_PROPEP(OF, CP)
        pypropep_data.append({
            'OF': f"{OF:.1f}", #just to make the indexing ok
            'CP': CP,
            'chamber_temp': chamber_temp,
            'molar_weight': molar_weight,
            'heat_ratio': heat_ratio
        })

ppp_df = pd.DataFrame(pypropep_data)
ppp_df.set_index(['OF', 'CP'], inplace=True)
#to get values, do ppp_df.loc[('1.0', 1200)] with 1st index as a string!!!

#ppp_df.to_json("backend\\PROPEP_lookup_table\\pypropep_lookup_table.json", orient="table", indent = 2)
#ppp_df_loaded = pd.read_json("backend\\PROPEP_lookup_table\\pypropep_lookup_table.json", orient="table")