
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from math import pi


# navigate to folder where this script is stored
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


fname = 'HOTFIRE2.7.txt'

#current testsite setup 2023-02-17
column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                'purge_pressure_sw', 'tank_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw','ven_temp_sw']


#for old data analysis 
#column_names = ['time_ms','raw_fill_pt', 'raw_tank_pt', 'raw_purge_pt', 'raw_run_pt','thrust_V', 
          #'raw_cc_pt','raw_tank_mt', 'raw_run_tt', 'raw_tank_tt', 'cjc', 'fill_pressure_sw', 'tank_pressure_sw', 
          #'purge_pressure_sw', 'run_pressure_sw','thrust_sw', 'cc_pressure_sw','tank_mt', 'run_tt', 'tank_tt']


# read file as dataframe
df = pd.read_csv(fname, sep='\t', names=column_names, index_col=False, dtype=np.float64, skip_blank_lines=True)



# start and end of actual hotfire data amid the huge file

hf_samples_start = 60000
hf_samples_end = 65000

#start and end of filling

fill_start = 0
fill_end = 950000000


# extract columns
time = df.loc[:, 'time_ms'] / 1000
time_zeroed = time - time[hf_samples_start]
thrust_raw = df.loc[:, 'thrust_V']
tank_temp = df.loc[:, 'tank_temp_sw']
run_temp = df.loc[:, 'run_temp_sw']
thrust = df.loc[:,'thrust_sw']
tank_mass = df.loc[:,'tank_mass_sw']

fill_pressure = df.loc[:,'fill_pressure_sw']
run_pressure = df.loc[:, 'run_pressure_sw'] 
cc_pressure = df.loc[:,'cc_pressure_sw']
tank_pressure = df.loc[:,'tank_pressure_sw']

#for rocket cold flow

tank_pressure_top = []
for i in thrust_raw:
    p = ((i-thrust_raw[hf_samples_start])*-47 + 500)
    tank_pressure_top.append(p)
    i+=1



gain = tank_pressure[hf_samples_start]/cc_pressure[hf_samples_start]
print(gain)

cc_pressure_new = []
for i in cc_pressure:
    p = ((i-cc_pressure[hf_samples_start])*-0.19 + 500)
    cc_pressure_new.append(p)
    i+=1

# k = 1.303 #specific heat ratio
# nozzle_radius = 0.0233 #nozzle throat radius
# A = pi*(nozzle_radius)**2 #nozzle throat area

# cc_pressure = []
# for i in thrust:
#     #calculate pressure in 
#     r = (((2*k**2)/(k-1))*((2)/(k+1))**((k+1)/(k-1)))**(1/2)
#     pressure_pascals = i/(A*r)
#     pressure_psi = pressure_pascals/6895
#     cc_pressure.append(pressure_psi)
#     i+=1

#get delta P plots
delta_P_tank = tank_pressure - run_pressure
delta_P_injector = tank_pressure_top - tank_pressure

#plot important pressure plots for hotfire

plt.figure()
#plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(tank_pressure[hf_samples_start : hf_samples_end]), label = 'Tank Pressure [psi]')
#plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(run_pressure[hf_samples_start : hf_samples_end]), 'r', label = 'Run Pressure [psi]')
#plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(delta_P_tank[hf_samples_start : hf_samples_end]), 'g--', label = 'Delta P [psi]')
plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(thrust_raw[hf_samples_start : hf_samples_end]), 'g--', label = 'Delta P [psi]')
plt.xlabel('Time (s)')  # Add an x-label to the axes.
plt.ylabel('Pressure [psi]')  # Add a y-label to the axes.
plt.title('Rocket Cold Flow Important Pressure Plots (Tank Emptying)')  # Add a title to the axes.
plt.minorticks_on()
plt.grid(True, which = 'major', color = 'black', linestyle = '-')
plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
plt.legend()  # Add a legend
plt.show()

# plt.figure()
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(thrust[hf_samples_start : hf_samples_end]), label = 'Thrust [N]')
# plt.xlabel('Time (s)')  # Add an x-label to the axes.
# plt.ylabel('Thrust [N]')  # Add a y-label to the axes.
# plt.title('Hotfire Attempt 2.5 Thrust')  # Add a title to the axes.
# plt.minorticks_on()
# plt.grid(True, which = 'major', color = 'black', linestyle = '-')
# plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
# plt.legend()  # Add a legend
# plt.show()

# plt.figure()
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(tank_pressure_top[hf_samples_start : hf_samples_end]), label = 'Top of Tank Pressure [psi]')
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(cc_pressure_new[hf_samples_start : hf_samples_end]), 'r', label = 'Bottom of Tank Pressure 1 (cc pressure transducer) [psi]')
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(tank_pressure[hf_samples_start : hf_samples_end]), 'g-', label = 'Bottom of Tank Pressure 2 (tank pressure transducer) [psi]')
# plt.xlabel('Time (s)')  # Add an x-label to the axes.
# plt.ylabel('Pressure [psi]')  # Add a y-label to the axes.
# plt.title('Rocket Coldflow Tank Pressure')  # Add a title to the axes.
# plt.minorticks_on()
# plt.grid(True, which = 'major', color = 'black', linestyle = '-')
# plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
# plt.legend()  # Add a legend
# plt.show()

# plt.figure()
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(tank_temp[hf_samples_start : hf_samples_end]), label = 'Tank Temperature [C]')
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(run_temp[hf_samples_start : hf_samples_end]), 'r', label = 'Run Temperature [C]')
# plt.xlabel('Time (s)')  # Add an x-label to the axes.
# plt.ylabel('Temperature [C]')  # Add a y-label to the axes.
# plt.title('Hotfire Attempt 2.5 Temperature Plots')  # Add a title to the axes.
# plt.minorticks_on()
# plt.grid(True, which = 'major', color = 'black', linestyle = '-')
# plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
# plt.legend()  # Add a legend
# plt.show()

# plt.figure()
# plt.plot(time_zeroed[hf_samples_start : hf_samples_end], np.array(tank_mass[hf_samples_start : hf_samples_end]), label = 'Tank Mass [kg]')
# plt.xlabel('Time (s)')  # Add an x-label to the axes.
# plt.ylabel('Tank Mass [kg]')  # Add a y-label to the axes.
# plt.title('Hotfire Attempt 2.5 Tank Emptying Mass')  # Add a title to the axes.
# plt.minorticks_on()
# plt.grid(True, which = 'major', color = 'black', linestyle = '-')
# plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
# plt.legend()  # Add a legend
# plt.show()


# #plot important pressure plots during filling
# plt.figure()
# plt.plot(time_zeroed[fill_start : fill_end], np.array(fill_pressure[fill_start : fill_end]), label = 'Fill Pressure [psi]')
# plt.plot(time_zeroed[fill_start : fill_end], np.array(tank_pressure[fill_start : fill_end]), label = 'Tank Pressure [psi]')
# plt.xlabel('Time (s)')  # Add an x-label to the axes.
# plt.ylabel('Pressure [psi]')  # Add a y-label to the axes.
# plt.title('Hotfire Attempt 2.6 (Misfire) Important Pressure Plots During Tank Filling')  # Add a title to the axes.
# plt.minorticks_on()
# plt.grid(True, which = 'major', color = 'black', linestyle = '-')
# plt.grid(True, which ='minor', color='grey', linestyle='--', linewidth = 0.5)
# plt.legend()  # Add a legend
# plt.show()

# impulse = 0
# i = hf_samples_start

# while i < hf_samples_end:
#     impulse += (((thrust[i]))*(time[i+1] - time[i]))
#     i += 1 
# print("total impulse is:", impulse)

