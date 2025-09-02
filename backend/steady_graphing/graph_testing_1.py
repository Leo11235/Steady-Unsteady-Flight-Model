import plotly.graph_objects as go
import plotly.express as px
import json

file = "./output_files/simulation_output.json"
with open(file, 'r') as f:
    data = json.load(f)

flight_dict = data["flight dict"]

# 3D scatter plot example
fig = go.Figure(data=[go.Scatter3d(
    x=flight_dict['time'],
    y=flight_dict['altitude'], 
    z=flight_dict['velocity'],
    mode='markers',
    marker=dict(size=3, color=flight_dict['acceleration'], colorscale='Viridis')
)])

fig.update_layout(scene=dict(
    xaxis_title='Time (s)',
    yaxis_title='Altitude (m)',
    zaxis_title='Velocity (m/s)'
))
fig.show()