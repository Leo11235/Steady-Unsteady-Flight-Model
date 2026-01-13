from pathlib import Path
import pandas as pd 
import matplotlib.pyplot as plt

from src.backend import validation_utils
from src.backend.steady import steady_main

_ROOT_DIR = Path(__file__).resolve().parent

def steady_validation():
    cur_valid = validation_utils.valid_hotfires
    rh = validation_utils.registered_hotfires
    categories = ["Burntime", "Total Impulse", "Peak Thrust", "Average Thrust"]
    category_keys = ["burntime", "total_impulse", "peak_thrust", "avg_thrust"]
    formats = ["7.2f", "7.1f", "7.2f", "7.2f"]
    label_snippet = "|   Sim   —   Act   —   Err   "#.replace("—", " ")

    print("\n\n     \u001B[1m", end="")
    for s in categories:
        lc = len(label_snippet) - len(s)
        print(" " * int((lc + 1) / 2), s,  " " * int(lc / 2), sep = "", end="")
    print("\u001B[0m")
    print("    ", label_snippet * len(categories), "|", sep="")
    print("-" * 4, ("|" + "-" * (len(label_snippet) - 1)) * len(categories), "|", sep="")

    for i in range(len(rh)):
        if (not rh[i] in cur_valid): continue
        # if (rh[i] != "4.1"): continue
        hotfire_path = f"{_ROOT_DIR}/src/backend/static_data/hotfires/real_data/HOTFIRE{rh[i]}.jsonc" #hotfire outputs
        steady_sim_path = f"{_ROOT_DIR}/src/backend/static_data/hotfires/simulation_configs/Hotfire_{rh[i]}_steady.jsonc" #sim inputs

        hotfire_results = validation_utils.graph_all(hotfire_path, i)
        steady_input = steady_main.get_rocket_inputs(steady_sim_path)
        steady_input['chamber pressure'] = hotfire_results['avg_cc_pressure'] / 0.00014504 # psi to pascal
        steady_dic = steady_main.run_with_inputs(steady_input, silent=True)
        steady_results = validation_utils.graph_steady(steady_dic, i)

        # print(f"HOTFIRE:\n\tavg_cc_pressure = {hotfire_results['avg_cc_pressure']}\n\tburntime: {hotfire_results['burntime']} (HARDCODED)\n\tpeak_thrust = {hotfire_results['max_thrust']}\n\tavg_thrust = {hotfire_results['avg_thrust']}\n\ttotal_impulse = {hotfire_results['total_impulse']}")
        # print(f"Steady Simulation:\n\tburntime: {steady_results['burntime']}\n\tavg_thrust = {steady_results['avg_thrust']}\n\ttotal_impulse = {steady_results['total_impulse']}")
        # plt.show()

        print(f"{rh[i]} ", sep="", end="")
        for i in range(len(categories)):
            hotfire_val = hotfire_results[category_keys[i]]
            steady_val = steady_results[category_keys[i]]
            error = 100 * (steady_val / hotfire_val - 1)

            hs = format(hotfire_val, formats[i])
            ss = format(steady_val, formats[i])
            err = format(error, "6.2f")

            print(f"| {hs}   {ss}   {err}% ", end="")
        print("|")
    print("\n")


if __name__ == "__main__":
    steady_validation()