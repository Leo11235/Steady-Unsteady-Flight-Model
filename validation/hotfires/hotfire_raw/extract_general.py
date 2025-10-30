import os, json

#current testsite setup 2023-02-17
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'purge_pressure_sw', 'tank_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw','ven_temp_sw']

#from 4.1
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'tank_pressure_sw', 'purge_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw', 'ven_temp_sw']


def convert_file (path):
    dic = {}
    with open(path, 'r') as file:
        length = 0
        for line in file:
            len_split = int(len(str(line.strip()).split()))
            if len_split != 0:
                length = len_split
                break
        
        print(f"{path} has words {length}")
        for i in range(length):
            dic[f"col{i}"] = []
        for line in file:
            split_arr = str(line.strip()).split()
            if int(len(split_arr)) == 0: continue
            # print (split_arr)
            for i in range(length - 1):
                # print(f"{i}: {split_arr[i]}")
                dic[f"col{i}"].append(split_arr[i])
    with open(f"{path[:-4]}.jsonc", "w") as f:
        json.dump(dic, f, indent=4)


if __name__ == "__main__":
    arr = sorted(os.listdir("validation/hotfires/hotfire_raw/"))
    for file in arr:
        if (file.endswith(".txt")):
            convert_file("validation/hotfires/hotfire_raw/" + file)