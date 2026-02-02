import json
import matplotlib.pyplot as plt
import numpy as np
from ..steady.steady_main import steady_main

registered_hotfires = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]
valid_hotfires = list(["3.1", "3.4", "3.5", "4.1"])
start_arr =    [0, 0, 0, 0, 0, 22.75, 9, 5.63, 5]
burntime_arr = [0, 0, 0, 0, 0, 10.5, 7.2, 8, 10]

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
    # plt.title(path.split("/")[-1] + ' cc_pressure graph')

    avg_cc_pressure = 0
    for i in range(len(pressure) - 1): 
        if (sec[i] < start_arr[J]): continue
        elif (sec[i] > start_arr[J] + burntime_arr[J]): break
        avg_cc_pressure += pressure[i] * (sec[i + 1] - sec[i])
    avg_cc_pressure /= burntime_arr[J]

    dic['avg_cc_pressure'] = avg_cc_pressure


def graph_thrust (path, sec, dict_k, dic, J):
    burntime = burntime_arr[J]
    starttime = start_arr[J]

    plt.figure('thrust')
    plt.xlabel('seconds')
    plt.ylabel('thrust')
    # plt.plot([sec[0], sec[-1]], [0, 0], color="black", linewidth=1) # zero axis

    total_impulse_uncorrected = 0
    for i in range(len(dict_k) - 1): total_impulse_uncorrected += dict_k[i] * (sec[i + 1] - sec[i])

    last_avg = avg_last_vals(dict_k, 50)
    first_index = getFirstLocAtVal(dict_k, last_avg)
    for i in range(len(dict_k)):
        if (i > first_index and dict_k[i] > 0.7 * last_avg): dict_k[i] -= last_avg
    
    #plt.plot(sec, dict_k) #adjusted thrust so it starts and ends ~0N
    #plot in three segments
    temp_thrusts = [[], [], []]
    temp_seconds = [[], [], []]
    for i in range(len(sec)):
        if (sec[i] < starttime): sel = 0
        elif (sec[i] < starttime + burntime): sel = 1
        else: sel = 2
        temp_thrusts[sel].append(dict_k[i])
        temp_seconds[sel].append(sec[i])
    plt.plot(temp_seconds[0], temp_thrusts[0], linewidth=0.5, color="gray")
    plt.plot(temp_seconds[2], temp_thrusts[2], linewidth=0.5, color="gray")
    plt.plot(temp_seconds[1], temp_thrusts[1], linewidth=1.2, color="orange")

    plt.title(path.split("/")[-1] + ' thrust graph')
    diff = (max(dict_k)-min(dict_k))
    if (diff == 0): diff = 0.1
    plt.yticks(np.arange(min(dict_k), max(dict_k), diff/10)) 

    tmax = 0
    for thrustpoint in dict_k:
        if (thrustpoint > tmax): tmax = thrustpoint
    total_impulse = 0
    for i in range(len(dict_k) - 1): total_impulse += dict_k[i] * (sec[i + 1] - sec[i])

    dic['peak_thrust'] = tmax
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
    pwlx = [0, start_time - margin, start_time, start_time + burntime, start_time + burntime+margin]
    pwly = [0, 0, thrust, thrust, 0]
    plt.plot(pwlx, pwly, color="darkgreen", linestyle="dashed", linewidth=1)

    total_impulse = steady_dic[1]['total impulse']
    # print(f"  steady_sim\n\ttotal_impulse = {total_impulse}\n\tthrust = {thrust}")

    dic = {}
    dic['avg_thrust'] = thrust
    dic['peak_thrust'] = thrust #to simplify print code, since these are the same in steady.
    dic['burntime'] = burntime
    dic['total_impulse'] = total_impulse
    return dic

#not tested yet, hopefully it works
def graph_unsteady(unsteady_dic):
    time = unsteady_dic["time"]
    F_x = unsteady_dic["F_x"]
    F_y = unsteady_dic["F_y"]
    total_thrust = np.sqrt(F_x**2 + F_y**2)

    threshold = 0.01 * np.max(total_thrust) #thrust cutoff threshold
    burn_idx = np.where(total_thrust >= threshold)[0]

    burntime_start = time[burn_idx[0]]
    burntime_end = time[burn_idx[-1]]
    burntime = burntime_end - burntime_start

    total_impulse = np.trapz(total_thrust[burn_idx], time[burn_idx])

    plt.figure("thrust")
    plt.plot(time, total_thrust, color="darkgreen", linestyle="dashed", linewidth=1)

    result_dic = {'avg_thrust': np.mean(total_thrust[burn_idx]),'peak_thrust': np.max(total_thrust),'burntime': burntime,'total_impulse': total_impulse}
    return result_dic


# if __name__ == "__main__":
#     cur_valid = list(["3.1", "3.2", "3.4", "3.5", "4.1"]);
#     for i in range(len(arr)):
#         if (not arr[i] in cur_valid): continue
#         if (arr[i] != "4.1"): continue
#         p = f"validation/hotfires/hotfire_processed/HOTFIRE{arr[i]}.jsonc"

#         hotfire_results = graph_all(p, i)
#         steady_input = get_rocket_inputs(f"steady_input_files/Hotfire_{arr[i]}.jsonc")
#         steady_input['chamber pressure'] = hotfire_results['avg_cc_pressure'] / 0.00014504 # psi to pascal
#         steady_dic = run_with_inputs(steady_input)
#         steady_results = graph_steady(steady_dic, i)

#         print(f"HOTFIRE:\n\tavg_cc_pressure = {hotfire_results['avg_cc_pressure']}\n\tburntime: {hotfire_results['burntime']}\n\tpeak_thrust = {hotfire_results['max_thrust']}\n\tavg_thrust = {hotfire_results['avg_thrust']}\n\ttotal_impulse = {hotfire_results['total_impulse']}")
#         print(f"Steady Simulation:\n\tburntime: {steady_results['burntime']}\n\tavg_thrust = {steady_results['avg_thrust']}\n\ttotal_impulse = {steady_results['total_impulse']}")
#         plt.show()