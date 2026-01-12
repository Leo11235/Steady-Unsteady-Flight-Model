'''
This file functions as a tester for venvs. Your output should look like:

Running fuel mass convergence simulation
    Loop 1 of simulation, apogee achieved: 14910.5 meters
    Loop 2 of simulation, apogee achieved: 8103.0 meters
    Loop 3 of simulation, apogee achieved: 11720.4 meters
    Loop 4 of simulation, apogee achieved: 13368.5 meters
    Loop 5 of simulation, apogee achieved: 14153.5 meters
    Loop 6 of simulation, apogee achieved: 13764.4 meters
    Loop 7 of simulation, apogee achieved: 13567.2 meters
    Loop 8 of simulation, apogee achieved: 13666.0 meters
    LOOP 9 - FINAL APOGEE: 13715.2 meters (0.81256 meters off from target)
    
'''

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import backend.steady_main

if __name__ == '__main__':
    file_path = "./steady_input_files/Esteban's_Ancalagon.jsonc"
    output_data = backend.steady_main.main(file_path)
