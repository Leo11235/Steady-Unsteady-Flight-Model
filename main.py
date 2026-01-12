from pathlib import Path
import pandas as pd 
import matplotlib.pyplot as plt

from src.backend import validation_utils
from src.backend.steady import steady_main

_ROOT_DIR = Path(__file__).resolve().parent

def steady_validation():
    cur_valid = validation_utils.valid_hotfires
    rh = validation_utils.registered_hotfires
    for i in range(len(rh)):
        if (not rh[i] in cur_valid): continue
        # if (rh[i] != "4.1"): continue
        hotfire_path = f"{_ROOT_DIR}/src/backend/static_data/hotfires/real_data/HOTFIRE{rh[i]}.jsonc" #hotfire outputs
        steady_sim_path = f"{_ROOT_DIR}/src/backend/static_data/hotfires/simulation_configs/Hotfire_{rh[i]}_steady.jsonc" #sim inputs

        hotfire_results = validation_utils.graph_all(hotfire_path, i)
        steady_input = steady_main.get_rocket_inputs(steady_sim_path)
        steady_input['chamber pressure'] = hotfire_results['avg_cc_pressure'] / 0.00014504 # psi to pascal
        steady_dic = steady_main.run_with_inputs(steady_input)
        steady_results = validation_utils.graph_steady(steady_dic, i)

        print(f"HOTFIRE:\n\tavg_cc_pressure = {hotfire_results['avg_cc_pressure']}\n\tburntime: {hotfire_results['burntime']}\n\tpeak_thrust = {hotfire_results['max_thrust']}\n\tavg_thrust = {hotfire_results['avg_thrust']}\n\ttotal_impulse = {hotfire_results['total_impulse']}")
        print(f"Steady Simulation:\n\tburntime: {steady_results['burntime']}\n\tavg_thrust = {steady_results['avg_thrust']}\n\ttotal_impulse = {steady_results['total_impulse']}")
        plt.show()

if __name__ == "__main__":
    steady_validation()