import json
import matplotlib.pyplot as plt
import numpy as np
from backend.steady_main import get_rocket_inputs, run_with_inputs

arr = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]
start_arr =    [0, 0, 0, 0, 0, 22.75, 9, 5.63, 5]
burntime_arr = [0, 0, 0, 0, 0, 19, 0, 0, 20]

def dic_of (path):
    with open (path) as f:
        return json.load(f)
    
def avg_last_vals(list, length):
    val = 0
    for i in range(length): val += list[-i]
    return val / length

def getFirstLocAtVal (list, value):
    for i in range(len(list)): 
        if (list[i] >= value): 
            return i
    raise Exception("no value >= " + value)

def graph_all (path, J):
    dic_data = dic_of(path)
    sec = dic_data['seconds']
    offset = sec[0]
    for i in range(len(sec)): sec[i] -= offset

    dic = {} #tbr

    for k in dic_data.keys():
        if (k == "seconds"): continue
        if (len(sec) != len(dic_data[k])):
            print(f"{path} + key {k} has length {len(dic_data[k])}")
            continue

        if (k == "thrust"): graph_thrust(path, sec, dic_data[k], dic, J)
        elif (k == "cc_pressure"): graph_pressure(path, sec, dic_data[k], dic, J)
    return dic

def graph_pressure (path, sec, pressure, dic, J):
    # plt.figure('cc_pressure')
    # plt.xlabel('seconds')
    # plt.ylabel('cc_pressure')
    # plt.plot(sec, pressure) #adjusted thrust so it starts and ends ~0N
    # plt.title(path + ' cc_pressure graph')

    avg_cc_pressure = 0
    for i in range(len(pressure) - 1): 
        if (sec[i] < start_arr[J]): continue
        elif (sec[i] > start_arr[J] + burntime_arr[J]): break
        avg_cc_pressure += pressure[i] * (sec[i + 1] - sec[i])
    avg_cc_pressure /= burntime_arr[J]

    dic['avg_cc_pressure'] = avg_cc_pressure


def graph_thrust (path, sec, dict_k, dic, J):
    plt.figure('thrust')
    plt.xlabel('seconds')
    plt.ylabel('thrust')
    plt.plot([sec[0], sec[-1]], [0, 0])

    total_impulse_uncorrected = 0
    for i in range(len(dict_k) - 1): total_impulse_uncorrected += dict_k[i] * (sec[i + 1] - sec[i])

    last_avg = avg_last_vals(dict_k, 50)
    first_index = getFirstLocAtVal(dict_k, last_avg)
    for i in range(len(dict_k)):
        if (i > first_index and dict_k[i] > 0.7 * last_avg): dict_k[i] -= last_avg
    
    plt.plot(sec, dict_k) #adjusted thrust so it starts and ends ~0N

    plt.title(path + ' thrust graph')
    diff = (max(dict_k)-min(dict_k))
    if (diff == 0): diff = 0.1
    plt.yticks(np.arange(min(dict_k), max(dict_k), diff/10)) 

    tmax = 0
    for thrustpoint in dict_k:
        if (thrustpoint > tmax): tmax = thrustpoint
    total_impulse = 0
    for i in range(len(dict_k) - 1): total_impulse += dict_k[i] * (sec[i + 1] - sec[i])

    burntime = burntime_arr[J]

    dic['max_thrust'] = tmax
    # dic['avg_cc_pressure'] = 
    dic['burntime'] = burntime
    dic['avg_thrust'] = total_impulse / burntime
    dic['total_impulse'] = total_impulse
    dic['total_impulse_uncorrected'] = total_impulse_uncorrected

def graph_steady (steady_dic, i):
    start_time = start_arr[i]
    thrust = steady_dic[1]['thrust']
    burntime = steady_dic[1]['burntime']
    margin = 0.001
    plt.figure("thrust")
    plt.plot([start_time - margin, start_time, start_time + burntime, start_time + burntime+margin], [0, thrust, thrust, 0])

    total_impulse = steady_dic[1]['total impulse']
    # print(f"  steady_sim\n\ttotal_impulse = {total_impulse}\n\tthrust = {thrust}")

    dic = {}
    dic['avg_thrust'] = thrust
    dic['burntime'] = burntime
    dic['total_impulse'] = total_impulse
    return dic


if __name__ == "__main__":
    cur_valid = list(["3.1", "3.2", "3.4", "3.5", "4.1"]);
    for i in range(len(arr)):
        if (not arr[i] in cur_valid): continue
        if (arr[i] != "4.1"): continue
        p = f"validation/hotfires/hotfire_processed/HOTFIRE{arr[i]}.jsonc"

        hotfire_results = graph_all(p, i)
        steady_input = get_rocket_inputs(f"steady_input_files/Hotfire_{arr[i]}.jsonc")
        steady_input['chamber pressure'] = hotfire_results['avg_cc_pressure'] / 0.00014504 # psi to pascal
        steady_dic = run_with_inputs(steady_input)
        steady_results = graph_steady(steady_dic, i)

        print(f"HOTFIRE:\n\tavg_cc_pressure = {hotfire_results['avg_cc_pressure']}\n\tburntime: {hotfire_results['burntime']}\n\tpeak_thrust = {hotfire_results['max_thrust']}\n\tavg_thrust = {hotfire_results['avg_thrust']}\n\ttotal_impulse = {hotfire_results['total_impulse']}")
        print(f"Steady Simulation:\n\tburntime: {steady_results['burntime']}\n\tavg_thrust = {steady_results['avg_thrust']}\n\ttotal_impulse = {steady_results['total_impulse']}")
        plt.show()


