import json
import os

processed_hotfires = ["hotfire2.7"]



if __name__ == '__main__':
    chosen_hotfire = processed_hotfires[0]
    # run a steady sim on given input parameters, write to given output.
    os.system(f"cd ..; python3 backend/steady_main.py validation/hotfires/hotfire_processed/{chosen_hotfire}/{chosen_hotfire}_params.jsonc validation/steady/{chosen_hotfire}.jsonc")
    steady_output = json.load(f"steady/{chosen_hotfire}.jsonc")

    # for steady state, we want to contrast, total thrust produced/avg thrust, and specific impulse, burn time, two masses, temp, nozzle
    # for unsteady we will be contrasting all factors
    pass