import json, debug_output_printing
import steady_backend.graphing as graphing

file_name = "Ancalagon 2 -- (CP[100, 1000], Lf[7.87, 78.74]) -- (30-01-2025, 15-01-35).json"
#file_name = "Ancalagon 2 -- (CP[300.00066026400003, 700.001540616]) -- (11-02-2025, 17-37-28).json"
file_name = "Ancalagon 2 -- (Mo[1, 8], CP[100, 1000]) -- (30-01-2025, 17-46-26).json"
file_name = "oldformat2.json"

file_path = f"./static/saved_graphs/{file_name}"

if __name__ == "__main__":
    #file_path = f"./static/saved_graphs/{file_name}" + (".json" if not file_name.endswith(".json") else "")
    with open(file_path, 'r') as json_file:
        parametric_data = json.load(json_file)
    graphing.graph_param(parametric_data)
    
    
    
    