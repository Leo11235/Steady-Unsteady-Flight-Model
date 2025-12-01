import json
from backend.steady_main import main, print_dict

if __name__ == '__main__':
    # get file and run simulation
    HOTFIRE = "3.5"
    file = f"./steady_validation_input_files/Hotfire_{HOTFIRE}_steady.jsonc"
    rocket_inputs, rocket_parameters = main(file)
    
    print(f"Thrust: {rocket_parameters.get("thrust", ":(")}")
    print(f"Total impulse: {rocket_parameters.get("total impulse")}")
    print(f"Fuel mass: {rocket_parameters.get("fuel mass")}")
    print(f"Burntime: {rocket_parameters.get("burntime")}")
    
    
    # # estimate 'file size' of data
    # json_str = json.dumps(data)
    # size_bytes = len(json_str.encode('utf-8'))
    # size_kb = size_bytes / 1024
    # size_mb = size_kb / 1024
    # print(f"Estimated JSON size: {size_bytes} bytes ({size_kb:.2f} KB / {size_mb:.2f} MB)")

