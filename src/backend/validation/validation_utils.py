import json
import matplotlib.pyplot as plt
import numpy as np
from ..steady.steady_main import steady_main

registered_hotfires = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]

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

def graph_hotfire_data (path):
    dic = dic_of(path)
    sec = dic['seconds']
    offset = sec[0]
    for i in range(len(sec)):
        sec[i] -= offset
    i = 0
    for k in dic.keys():
        if (k == "seconds"): continue
        if (len(sec) != len (dic[k])):
            print(f"----{path} + key {k} has length {len(dic[k])}")
            continue
        i = i + 1
        plt.figure(i - 1)
        plt.xlabel('seconds')
        plt.ylabel(k)
        # plt.plot(sec, dic[k]) # data thrust
        plt.plot([sec[0], sec[-1]], [0, 0])

        # total_impulse = 0
        # for i in range(len(dic[k]) - 1): total_impulse += dic[k][i] * (sec[i + 1] - sec[i])
        # print(f"{path}: total impulse (uncorrected) {total_impulse}")

        last_avg = avg_last_vals(dic[k], 50)
        first_index = getFirstLocAtVal(dic[k], last_avg)
        for i in range(len(dic[k])):
            if (i > first_index and dic[k][i] > 0.7 * last_avg): dic[k][i] -= last_avg
        
        plt.plot(sec, dic[k]) #adjusted thrust so it starts and ends ~0N

        total_impulse = 0
        for i in range(len(dic[k]) - 1): total_impulse += dic[k][i] * (sec[i + 1] - sec[i])
        hotfire_name = path.split("/")[-1]
        print(f"----{hotfire_name}: total impulse {total_impulse}")

        plt.title(hotfire_name + ' graph ' + k)
        diff = (max(dic[k])-min(dic[k]))
        if (diff == 0): diff = 0.1
        plt.yticks(np.arange(min(dic[k]), max(dic[k]), diff/10)) 

def graph_steady (steady_dic, start_time):
    # print(steady_dic)
    thrust = steady_dic[1]['thrust']
    burntime = steady_dic[1]['burntime']
    margin = 0.001
    plt.plot([start_time - margin, start_time, start_time + burntime, start_time + burntime+margin], [0, thrust, thrust, 0])

    total_impulse = steady_dic[1]['total impulse']
    print(f"----steady_sim total_impulse = {total_impulse}")


# if __name__ == "__main__":
#     for n in arr:
#         if (n != "4.1"): continue
#         p = f"validation/hotfires/hotfire_processed/HOTFIRE{n}.jsonc"
#         dic = dic_of(p)
#         graph_all(p)
#         graph_steady(steady_main(f"steady_validation_input_files/Hotfire_{n}_steady.jsonc"), 5)
#         plt.show()
#         break
