import os, json
import matplotlib.pyplot as plt
import numpy as np

#current testsite setup 2023-02-17
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'purge_pressure_sw', 'tank_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw','ven_temp_sw']

#from 4.1
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'tank_pressure_sw', 'purge_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw', 'ven_temp_sw']

references = [
    # name, start_t, end_t, thrust_col
    # ["2.1", -1, -1],
    ["2.2", 7380, 7400, "col16"],
    ["2.3", 4775, 4850, "col17"],
    ["2.4", 3600, 3820, "col17"],
    # ["2.5", -1, -1],
    ["2.6", 4300, 4375, "col17"], #col17 ? why so low
    ["2.7", 3450, 3525, "col17"],
    # ["2.8", -1, -1],
    ["3.1", 4675, 4735, "col17"],
    # ["3.2", 4570, 4610],
    # ["3.3", 4430, 4470],
    ["3.4", 4330, 4360, "col17"],
    ["3.5", 2200, 2240, "col17"],
    ["4.1", 6625, 6675, "col17"],# off by 1000
]


def convert_file (path, t_start, t_end):
    dic = {}
    dic_shortened = {}
    with open(path, 'r') as file:
        length = 0
        for line in file:
            len_split = int(len(str(line.strip()).split()))
            if len_split != 0:
                length = len_split
                break
        
        print(f"{path} has words {length}")
        colname = {}
        for i in range(length):
            colname[i] = f"col{i}"

            if (i == 0): colname[i] = "seconds";

            dic[colname[i]] = []
            dic_shortened[colname[i]] = []

        for line in file:
            split_arr = str(line.strip()).split()
            if int(len(split_arr)) == 0: continue
            # print (split_arr)
            for i in range(length - 1):
                # print(f"{i}: {split_arr[i]}")
                dic[colname[i]].append(float(split_arr[i]))

    seconds_offset = float(dic['seconds'][0]) / 1000
    for i in range(len(dic['seconds'])): 
        tval = float(dic['seconds'][i]) / 1000 - seconds_offset
        dic['seconds'][i] = tval
        if (t_start == -1 or tval > t_start and tval < t_end): 
            for X in range(length - 1): 
                dic_shortened[colname[X]].append(dic[colname[X]][i])
    
    with open(f"{path[:-4]}.jsonc", "w") as f:
        json.dump(dic_shortened, f, indent=4)

    c1 = dic_shortened['col1']
    plt.plot(dic_shortened['seconds'], c1)
    plt.xlabel('seconds')
    plt.ylabel('col1')
    plt.title(path + ' graph')
    plt.yticks(np.arange(min(c1), max(c1), 0.1))
    plt.show()

def print_all_graphs (path):
    dic = {}
    with open (path) as f:
        dic = json.load(f)

#from hotfire 4.1 
    # dic['thrust_sw'] *= 1000
    # dic['cc_pressure_sw'] = 0.0051 * dic['cc_pressure_V'] + 0.2499
    # dic[f"col{column_names.index('thrust_sw')}"] = [i * 1000 for i in dic[f"col{column_names.index('thrust_sw')}"]]
    # dic[f"col{column_names.index('cc_pressure_sw')}"] = [0.0051 * i + 0.2499 for i in dic[f"col{column_names.index('cc_pressure_V')}"]]

    i = 0
    sec = dic['seconds']
    for k in dic.keys():
        if (k == "seconds"): continue
        if (len(sec) != len (dic[k])):
            print(f"{path} + key {k} has length {len(dic[k])}")
            continue
        i = i + 1
        if (i < 16): continue
        # if (column_names[i].__contains__("_V")): continue
        plt.figure(i - 1)
        plt.xlabel('seconds')
        plt.ylabel(k)
        plt.plot(sec, dic[k])
        plt.title(path + ' graph ' + k + " " + column_names[i] + "???")
        diff = (max(dic[k])-min(dic[k]))
        if (diff == 0): diff = 0.1
        plt.yticks(np.arange(min(dic[k]), max(dic[k]), diff/10))
    plt.show()

def read_write_raw_processed_jsonc (path, reference_list):
    dic = {}
    with open (path) as f:
        dic = json.load(f)
    
    output_dic = {}
    output_dic['seconds'] = dic['seconds']
    output_dic['thrust'] = dic[reference_list[3]]
    if (reference_list[0] == "4.1"): output_dic['thrust'] = [1000 * i for i in output_dic['thrust']]

    with open(f"validation/hotfires/hotfire_processed/HOTFIRE{reference_list[0]}.jsonc", "w") as f:
        json.dump(output_dic, f, indent=4)

def graph_all (path):
    dic = {}
    with open (path) as f:
        dic = json.load(f)
    sec = dic['seconds']
    i = 0
    for k in dic.keys():
        if (k == "seconds"): continue
        if (len(sec) != len (dic[k])):
            print(f"{path} + key {k} has length {len(dic[k])}")
            continue
        i = i + 1
        plt.figure(i - 1)
        plt.xlabel('seconds')
        plt.ylabel(k)
        plt.plot(sec, dic[k])
        plt.title(path + ' graph ' + k)
        diff = (max(dic[k])-min(dic[k]))
        if (diff == 0): diff = 0.1
        plt.yticks(np.arange(min(dic[k]), max(dic[k]), diff/10))    
    plt.show()


if __name__ == "__main__":
    arr = reversed(references)
    for list in arr:
        # if (list[0] != "4.1"): continue
        # convert_file("validation/hotfires/hotfire_raw/HOTFIRE" + list[0] + ".txt", list[1], list[2])
        # print_all_graphs("validation/hotfires/hotfire_raw/HOTFIRE" + list[0] + ".jsonc")
        read_write_raw_processed_jsonc("validation/hotfires/hotfire_raw/HOTFIRE" + list[0] + ".jsonc", list)
        graph_all("validation/hotfires/hotfire_processed/HOTFIRE" + list[0] + ".jsonc")
        