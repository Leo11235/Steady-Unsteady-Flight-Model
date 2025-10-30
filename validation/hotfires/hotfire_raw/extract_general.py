import os

#current testsite setup 2023-02-17
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'purge_pressure_sw', 'tank_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw','ven_temp_sw']

#from 4.1
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'tank_pressure_sw', 'purge_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw', 'ven_temp_sw']


def convert_file (path):
    min = 100
    max = 0
    with open(path, 'r') as file:
        for line in file:
            stripped = str(line.strip())
            len_split = int(len(stripped.split()))
            if len_split != 0:
                if (len_split > max): max = len_split
                elif (len_split < min): min = len_split
        print(f"{path} has max words {max} and min words {min}")


if __name__ == "__main__":
    arr = sorted(os.listdir("validation/hotfires/hotfire_raw/"))
    for file in arr:
        if (file.endswith(".txt")):
            convert_file("validation/hotfires/hotfire_raw/" + file)