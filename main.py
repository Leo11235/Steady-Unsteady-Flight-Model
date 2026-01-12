from pathlib import Path
import pandas as pd 
import matplotlib.pyplot as plt

from src.backend import steady_main, unsteady_main, validation_utils

_ROOT_DIR = Path(__file__).resolve().parent

def steady_validation():
    for n in validation_utils.registered_hotfires:
        if (n != "4.1"): continue
        p = f"{_ROOT_DIR}/src/backend/static_data/hotfires/real_data/HOTFIRE{n}.jsonc"
        validation_utils.graph_hotfire_data(p)
        steady_output = steady_main(f"{_ROOT_DIR}/src/backend/static_data/hotfires/simulation_configs/Hotfire_{n}_steady.jsonc")
        validation_utils.graph_steady(steady_output, 5)
        plt.show()
        break

if __name__ == "__main__":
    steady_validation()