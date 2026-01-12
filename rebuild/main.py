from pathlib import Path
import pandas as pd 
import matplotlib.pyplot as plt

from src.backend import steady_main, unsteady_main, validation


_ROOT_DIR = Path(__file__).resolve().parent
_SIM_CONFIGS_PATH = _ROOT_DIR / "user_data" / "simulation_configs"


# print(_SIM_CONFIGS_PATH)
# print(steady_main(_SIM_CONFIGS_PATH / "steady_Esteban's_Ancalagon.jsonc"))


def steady_validation():
    for n in validation.arr:
        if (n != "4.1"): continue
        #p = f"validation/hotfires/hotfire_processed/HOTFIRE{n}.jsonc"
        p = f"{_ROOT_DIR}/src/backend/static_data/hotfires/real_data/HOTFIRE{n}.jsonc"
        dic = validation.dic_of(p)
        validation.graph_all(p)
        #validation.graph_steady(steady_main(f"steady_validation_input_files/Hotfire_{n}_steady.jsonc"), 5)
        validation.graph_steady(steady_main(f"{_ROOT_DIR}/src/backend/static_data/hotfires/simulation_configs/Hotfire_{n}_steady.jsonc"), 5)
        plt.show()
        break

steady_validation()