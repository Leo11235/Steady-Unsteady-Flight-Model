import os
import math
import datetime
import numpy as np
import pandas as pd
from math import pi
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


######################################################## 1. DEFINE FUNCTIONS ########################################################

################## 1.1 Extract Data ##################

def extract_data(file_name):
    #current testsite setup 2023-02-17
    column_names = ['time_ms', 'run_pressure_V', 'fill_pressure_V', 'purge_pressure_V', 'tank_pressure_V', 'tank_mass_V', 
                    'thrust_V', 'cc_pressure_V', 'tank_temp_V', 'run_temp_V', 'vent_temp_V', 'garbage', 'run_pressure_sw', 'fill_pressure_sw', 
                    'purge_pressure_sw', 'tank_pressure_sw', 'tank_mass_sw', 'thrust_sw', 'cc_pressure_sw', 'run_temp_sw', 'tank_temp_sw','ven_temp_sw']

    # read file as dataframe
    df = pd.read_csv(file_name, sep='\t', names=column_names, index_col=False, dtype=np.float64, skip_blank_lines=True)

    time = (df['time_ms'] / 1000).tolist()                  # seconds
    p_fill = df['fill_pressure_sw'].tolist()                # psi
    p_tank = df['tank_pressure_sw'].tolist()                # psi
    p_run = df['run_pressure_sw'].tolist()                  # psi
    p_cc = df['cc_pressure_sw'].tolist()                    # psi
    tank_temp = (df['tank_temp_sw']).tolist()               # celcius
    run_temp = (df['run_temp_sw']).tolist()                 # celcius
    tank_mass = (df['tank_mass_sw']).tolist()               # ???
    thrust = (df['thrust_sw']).tolist()                     # N
    mass = (df['tank_mass_sw']).tolist()                    # ???

    
    # Filling Time Range - To be automated
    start_time_fill = 5850
    end_time_fill = 7800
    time_range_fill = [t for t in time if start_time_fill <= t <= end_time_fill]
    time_zeroed_fill = [t - start_time_fill for t in time_range_fill]
    
    # Hotfire time range - To be automated
    start_time = 7800
    end_time = 7825  
    time_range = [t for t in time if start_time <= t <= end_time]
    time_zeroed = [t - start_time for t in time_range]
    
    # Full time range - To be automated
    start_time_full = 7800
    end_time_full = 7825
    time_range_full = [t for t in time if start_time_full <= t <= end_time_full]
    time_zeroed_full = [t - start_time_full for t in time_range_full]
    
    
    # combustion chamber processing
    p_cc_new = []
    for i in p_cc:
        p = ((i-p_cc[start_time])*0.19) #scaling of 0.19?
        p_cc_new.append(p)
        i+=1
    
    # get data ranges in time frame
    p_fill_range = [p_f for t, p_f in zip(time, p_fill) if start_time_fill <= t <= end_time_fill]
    p_tank_fill_range = [p_t_f for t, p_t_f in zip(time, p_tank) if start_time_fill <= t <= end_time_fill]
    p_tank_range = [p_t for t, p_t in zip(time, p_tank) if start_time <= t <= end_time]
    p_run_range = [p_r for t, p_r in zip(time, p_run) if start_time <= t <= end_time]
    p_cc_new_range = [p_c for t, p_c in zip(time, p_cc_new) if start_time <= t <= end_time]
    tank_temp_range = [t_t for t, t_t in zip(time, tank_temp) if start_time <= t <= end_time]
    run_temp_range = [r_t for t, r_t in zip(time, run_temp) if start_time <= t <= end_time]
    thrust_range = [thr for t, thr in zip(time, thrust) if start_time <= t <= end_time]
    mass_range = [mass for t, mass in zip(time, mass) if start_time_full <= t <= end_time_full]
    
    # calculate pressure deltas
    delta_p_tank_fill = [p_tank - p_fill for p_tank, p_fill in zip(p_tank_fill_range, p_fill_range)]
    delta_p_tank_run = [p_tank - p_run for p_tank, p_run in zip(p_tank_range, p_run_range)]
    delta_p_run_cc = [p_run - p_cc_new for p_run, p_cc_new in zip(p_run_range, p_cc_new_range)]
    
    return time_zeroed, time_zeroed_fill, time_zeroed_full, p_tank_fill_range, p_fill_range, p_tank_range, p_run_range, p_cc_new_range, tank_temp_range, run_temp_range, delta_p_tank_run, delta_p_run_cc, delta_p_tank_fill, thrust_range, mass_range



################## 1.2 Tank Mass (Estimated) - Incomplete ##################

## variables
p_c = 7251          # kPa
T_c = 309.57        # Kelvin
rho_c = 452         # kg/m^3

# 4.1 constants
b1 = -6.71893
b2 = 1.35966
b3 = -1.3779
b4 = -4.051

# 4.2 constants
b1_l = 1.72328
b2_l = -0.83950
b3_l = 0.51060
b4_l = -0.10412

# 4.3 constants
b1_g = -1.00900
b2_g = -6.28792
b3_g = 7.50332
b4_g = -7.90463
b5_g = 0.629427

## tank dimensions
outer_radius = 0.09                                           # m
inner_radius = 0.08268                                        # m
cyl_height = 0.685                                            # m
V_total = 0.016                                               # m^3 
V_inner_cyl = np.pi * (inner_radius**2) * cyl_height          # m^3
V_inner_hem = V_total - V_inner_cyl                           # m^3

def equations(temp_kelvin):
    T_r = temp_kelvin/T_c
    
    # 4.1 Vapour Pressure
    p = p_c * math.exp((1/T_r) * ((b1 * (1-T_r)) + (b2 * (1-T_r)**(3/2)) + (b3 * (1-T_r)**(5/2)) + (b4 * (1-T_r)**5)))
    
    # 4.2 Density of Saturated Liquid
    p_l = rho_c * math.exp((b1_l * (1-T_r)**(1/3)) + (b2_l * (1-T_r)**(2/3)) + (b3_l * (1-T_r)) + (b4_l * (1-T_r)**(4/3)))
    
    # 4.3 Density of Saturated Gas
    p_g = rho_c * math.exp((b1_g * ((1/T_r)-1)**(1/3)) + (b2_g * ((1/T_r)-1)**(2/3)) + (b3_g * ((1/T_r)-1)) + (b4_g * ((1/T_r)-1)**(4/3)) + (b5_g * ((1/T_r)-1)**(5/3)))
    
    return p, p_l, p_g

## create dictionary of vapour pressure and temperature values (as solving for temp from pressure is too mathematically intensive)
temp_list = np.linspace(183, 309, num = 10000)
vap_list, liquid_list, gas_list = equation_lists(temp_list)
pressure_temp_dict = dict(zip(vap_list, temp_list))


def liquid_nitrous_mass(outside_temp, dip_tube_length):
    
    p, p_liq, p_gas = equations(outside_temp + 273.15)           # get pressure and density values given known temperature
    V_cyl_liq_ratio = 1 - (dip_tube_length / cyl_height)         # get percentage of volume in cylinder taken by liquid given ullage height
    V_hem_ratio = V_inner_hem/V_total                            # get percentage of volume taken by hemisphere
    V_liq_ratio = V_cyl_liq_ratio + V_hem_ratio                  # get volume of liquid in tank ratio
    V_liq = V_liq_ratio * V_total                                # get volume of liquid
    V_gas = V_total - V_liq                                      # get volume of gas
    m_liq = V_liq * p_liq                                        # get mass of liquid in tank
    m_gas = V_gas * p_gas                                        # get mass of gas in tank
    
    
    print('With a dip tube length of', dip_tube_length, 'm, and a daily average ambient temperature of', outside_temp, 'degrees Celcius,')
    print('the mass of liquid nitrous in the run tank is', m_liq, 'kg')
    print('the mass of gas nitrous in the run tank is', m_gas, 'kg')
    print('The total mass of nitrous in the tank is therefore', m_liq + m_gas, 'kg')



################## 1.3 Mass Flow Rate (Estimated) - Incomplete ##################

def mass_flow_rate(C_i, num_holes, hole_diameter, p_run, p_cc_new, run_temp):
    # define variables
    #W_o = 7.3074 * (10)**(-25)
    W_o = 0.044013                           # kg/mol (oxidizer_molar_weight)
    c_1 = 2.781
    c_2 = 0.27244
    c_3 = 309.57
    c_4 = 0.2882
    
    # injection area
    A_i = (np.pi * ((hole_diameter/1000)/2)**2)

    # saturated molar volume as a fn of temp
    v_l = []
    for temp in run_temp_range:
        temp_K = temp + 273.15
        v = ((c_2)**(1 + (1 - (temp_K/c_3))**c_4)) / c_1
        v_l.append(v * 1000)

    # mass flow_rate
    n_o = []
    for i in range(len(p_run_range)):
        n = C_i * num_holes * A_i * math.sqrt((2 * p_run_range[i] - p_cc_new_range[i]) / (W_o * v_l[i]))
        n_o.append(n)
    
    return n_o


################## 1.4 Generate Report ##################

def generate_report(template_filename, title, date, abstract_text, analysis_text, conclusions_text):
    # Read the LaTeX template file
    with open(template_filename, 'r') as f:
        template_content = f.read()

    # Define consistent data
    author = 'Propulsion Subteam - McGill Rocket Team'
    issues_table = ''

    # Replace placeholders in the template with the dynamic content
    report_content = template_content.replace('{{ title }}', title)
    report_content = report_content.replace('{{ author }}', author)
    report_content = report_content.replace('{{ date }}', date)
    report_content = report_content.replace('{{ abstract }}', abstract_text)
    report_content = report_content.replace('{{ analysis }}', analysis_text)
    report_content = report_content.replace('{{ issues_table }}', issues_table)
    report_content = report_content.replace('{{ conclusions }}', conclusions_text)

    # Write the modified content to a new .tex file
    with open('report.tex', 'w') as f:
        f.write(report_content)

    print(f"Report generated successfully. You can find it in 'report.tex'.")





######################################################## 2. Generate Plots ########################################################

################## 2.1 Extract Data ##################
file_name = 'HOTFIRE2.7.txt'
time_zeroed, time_zeroed_fill, time_zeroed_full, p_tank_fill_range, p_fill_range, p_tank_range, p_run_range, p_cc_new_range, tank_temp_range, run_temp_range, delta_p_tank_run, delta_p_run_cc, delta_p_tank_fill, thrust_range, mass_range = extract_data(file_name)


################## 2.2 Fill Pressure ##################
plt.plot(time_zeroed_fill, p_fill_range, label='Fill Pressure')
plt.plot(time_zeroed_fill, p_tank_fill_range, label='Tank Pressure')
plt.plot(time_zeroed_fill, delta_p_tank_fill, label='Delta_P')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Fill Pressure')
plt.grid(True)
plt.legend()
plt.savefig('tank_fill_pressure_drop.png')
plt.show()

################## 2.3 Temperature ##################
plt.plot(time_zeroed, tank_temp_range, label='Tank Temperature (C)')
plt.plot(time_zeroed, run_temp_range, label='Run Temperature (C)')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (C)')
plt.title('Fill Pressure')
plt.grid(True)
plt.legend()
plt.savefig('temp.png')
plt.show()

################## 2.4 Tank/Run Pressure ##################
plt.plot(time_zeroed, p_tank_range, label='Tank Pressure')
plt.plot(time_zeroed, p_run_range, label='Run Line Pressure')
plt.plot(time_zeroed, delta_p_tank_run, label='Delta_P')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Tank Emptying Pressure Drop')
plt.grid(True)
plt.legend()
plt.savefig('tank_run_pressure_drop.png')
plt.show()

################## 2.5 Run/CC Pressure ##################
plt.plot(time_zeroed, p_run_range, label='Run Line Pressure')
plt.plot(time_zeroed, p_cc_new_range, label='Combustion Chamber Pressure')
plt.plot(time_zeroed, delta_p_run_cc, label='Delta_P')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Injector Pressure Drop')
plt.grid(True)
plt.legend()
plt.savefig('run_cc_pressure_drop.png')
plt.show()

################## 2.6 Thrust ##################

plt.plot(time_zeroed, thrust_range, label='Thrust (N)')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (N)')
plt.title('Thrust')
plt.grid(True)
plt.legend()
plt.savefig('thrust.png')
plt.show()

################## 2.7 Mass (Measured) ##################
plt.plot(time_zeroed_full, mass_range, label='Mass ()')
plt.xlabel('Time (s)')
plt.ylabel('Mass ()')
plt.title('Tank Mass (Load Cell)')
plt.grid(True)
plt.legend()
plt.savefig('mass.png')
plt.show()

################## 2.8 Mass (Estimated) ##################
plt.plot(time_zeroed_fill, mass_range_estimate, label='Mass ()')
plt.xlabel('Time (s)')
plt.ylabel('Mass ()')
plt.title('Tank Mass (Estimated)')
plt.grid(True)
plt.legend()
plt.savefig('mass.png')
plt.show()


################## 2.9 Mass Flow Rate ##################
C_i = 0.65                                # discharge coefficient
num_holes = 22
hole_diameter = 2  
n_o = mass_flow_rate(C_i, num_holes, hole_diameter, p_run, p_cc_new, run_temp)

plt.plot(time_zeroed, n_o, label='Mass Flow Rate')
plt.xlabel('Time (s)')
plt.ylabel('Mass Flow Rate')
plt.title('Injector Mass Flow Rate')
plt.grid(True)
plt.legend()
plt.savefig('mass_flow_rate.png')
plt.show()


######################################################## 2. Generate Report ########################################################
title = 'HotFire 2.7'
date = '2023/03/11'
abstract_text = 'The purpose of this hotfire was to validate engine performance and validate that that the engine can be fired twice in one day. Initial run tank pressure of ... Peak pre-injector pressure of ... Peak combustion chamber pressure ... Estimated nitrous mass ... The peak thrust was ... with a total impulse of ... The mass flow rate was...'
analysis_text = 'Specific analysis of hotfire data to determine whether the specific goal was met'
conclusions_text = 'Here you briefly summarize your findings and next steps...'

# Generate the latex report
generate_report('hotfire_template.tex', title, date, abstract_text, analysis_text, conclusions_text)
