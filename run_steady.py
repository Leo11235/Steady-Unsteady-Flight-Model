import json
from backend.steady_main import main, print_dict

if __name__ == '__main__':
    # get file and run simulation
    # file = "./steady_validation_input_jsons/Joel's_steady_comparison_1.jsonc"
    # data = main(file)
    # print_dict(data["rocket parameters"])
    
    file = "./steady_validation_input_jsons/Hotfire_4.1_steady.jsonc"
    data = main(file)
    print_dict(data)
    
    # estimate 'file size' of data
    json_str = json.dumps(data)
    size_bytes = len(json_str.encode('utf-8'))
    size_kb = size_bytes / 1024
    size_mb = size_kb / 1024
    print(f"Estimated JSON size: {size_bytes} bytes ({size_kb:.2f} KB / {size_mb:.2f} MB)")

