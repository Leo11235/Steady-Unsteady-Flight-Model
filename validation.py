import json
import matplotlib.pyplot as plt
import numpy as np
from backend.steady_main import main, print_dict

arr = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]
start_arr = [0, 0, 0, 0, 0, 22.75, 9, 5.63, 5]

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

def graph_all (path):
    dic = dic_of(path)
    sec = dic['seconds']
    offset = sec[0]
    for i in range(len(sec)):
        sec[i] -= offset
    i = 0
    for k in dic.keys():
        if (k == "seconds"): continue
        if (len(sec) != len (dic[k])):
            print(f"{path} + key {k} has length {len(dic[k])}")
            continue
        i = i + 1
        plt.figure(i - 1)
        graph_thrust(path, sec, dic[k])

def graph_thrust (path, sec, dict_k):
    plt.xlabel('seconds')
    plt.ylabel('thrust')
    # plt.plot(sec, dict_k) # data thrust
    plt.plot([sec[0], sec[-1]], [0, 0])

    # total_impulse = 0
    # for i in range(len(dict_k) - 1): total_impulse += dict_k[i] * (sec[i + 1] - sec[i])
    # print(f"{path}: total impulse (uncorrected) {total_impulse}")

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
    print(f"{path}:\n\ttotal impulse {total_impulse} Ns\n\tpeak thrust = {tmax}N")

def graph_steady (steady_dic, start_time):
    thrust = steady_dic[1]['thrust']
    burntime = steady_dic[1]['burntime']
    margin = 0.001
    plt.plot([start_time - margin, start_time, start_time + burntime, start_time + burntime+margin], [0, thrust, thrust, 0])

    total_impulse = steady_dic[1]['total impulse']
    print(f"  steady_sim\n\ttotal_impulse = {total_impulse}\n\tthrust = {thrust}")


if __name__ == "__main__":
    for i in range(len(arr)):
        if (arr[i] != "3.5"): continue
        p = f"validation/hotfires/hotfire_processed/HOTFIRE{arr[i]}.jsonc"
        dic = dic_of(p)
        graph_all(p)
        graph_steady(main(f"steady_inputs_files/Hotfire_{arr[i]}.jsonc"), start_arr[i])
        plt.show()
        break
