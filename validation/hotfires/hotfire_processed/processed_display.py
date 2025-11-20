import os, json
import matplotlib.pyplot as plt
import numpy as np
import backend as bk

arr = ["2.2", "2.3", "2.4", "2.6", "2.7", "3.1", "3.4", "3.5", "4.1"]

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
        print(f"{path}: total impulse {total_impulse}")

        plt.title(path + ' graph ' + k)
        diff = (max(dic[k])-min(dic[k]))
        if (diff == 0): diff = 0.1
        plt.yticks(np.arange(min(dic[k]), max(dic[k]), diff/10))    

if __name__ == "__main__":
    for n in arr:
        p = f"validation/hotfires/hotfire_processed/HOTFIRE{n}.jsonc"
        dic = dic_of(p)
        # start, end = get_start_end_spike(2, dic['thrust'])
        # print(dic['seconds'][start], dic['seconds'][end])
        graph_all(p)
        steady = bk.steady_main.main("./steady_input_files/Esteban's_Ancalagon.jsonc")
        # plt.plot([dic['seconds'][0], dic['seconds'][len(dic['seconds']) - 1]], [steady[]])
        plt.show()
