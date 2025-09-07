import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
import mplcursors






def graph_param(parametric_data):
    #print(len(parametric_data["parametric inputs"]))
    plot_one_input_vs_Isp(parametric_data) if len(parametric_data["parametric inputs"]) == 5 else None
    plot_two_inputs_vs_Isp(parametric_data) if len(parametric_data["parametric inputs"]) == 9 else None
    #plot_three_inputs_vs_Isp(parametric_data) if len(parametric_data["parametric inputs"]) == 9 else None
    


# ------------------------------------------------------------
# -------------------- PLOT ROCKET ASCENT --------------------
# ------------------------------------------------------------

def graph_everything(timesteps_dict):
    """Plot position, velocity, acceleration, and forces in separate windows simultaneously."""
    plot_position(timesteps_dict, block=False)
    plot_velocity(timesteps_dict, block=False)
    plot_acceleration(timesteps_dict, block=False)
    plot_forces(timesteps_dict, block=False)
    plt.show()

def plot_position(timesteps_dict, block=False):
    """Plot position as a function of time."""
    time = timesteps_dict["time"]
    position = timesteps_dict["altitude"]
    plt.figure()
    plt.plot(time, position, label="Altitude", linewidth=2)
    plt.title("Rocket Altitude vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Position (m)")
    plt.grid(True)
    plt.legend()

def plot_velocity(timesteps_dict, block=False):
    """Plot velocity as a function of time."""
    time = timesteps_dict["time"]
    velocity = timesteps_dict["velocity"]
    plt.figure()
    plt.plot(time, velocity, label="Velocity", linewidth=2, color="orange")
    plt.title("Rocket Velocity vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.grid(True)
    plt.legend()

def plot_acceleration(timesteps_dict, block=False):
    """Plot acceleration as a function of time."""
    time = timesteps_dict["time"]
    acceleration = timesteps_dict["acceleration"]
    plt.figure()
    plt.plot(time, acceleration, label="Acceleration", linewidth=2, color="green")
    plt.title("Rocket Acceleration vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.grid(True)
    plt.legend()

def plot_forces(timesteps_dict, block=False):
    """Plot net force, drag force, and gravitational force as functions of time."""
    time = timesteps_dict["time"]
    net_force = timesteps_dict["net force"]
    drag_force = timesteps_dict["drag force"]
    gravitational_force = timesteps_dict["grav force"]
    thrust_force = timesteps_dict["thrust"]

    plt.figure()
    plt.plot(time, net_force, label="Net Force", linewidth=2, color="blue")
    plt.plot(time, thrust_force, label="Thrust", linewidth=2, color="orange")
    plt.plot(time, drag_force, label="Drag Force", linewidth=2, color="red")
    plt.plot(time, gravitational_force, label="Gravitational Force", linewidth=2, color="purple")
    plt.title("Rocket Forces vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.grid(True)
    plt.legend()

    
# -------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY GRAPH: 1 INPUT --------------------
# -------------------------------------------------------------------------

def plot_one_input_vs_Isp(parametric_data):
    # unpack data
    graph_title, inputs, outputs = unpack_data(parametric_data)
    # get x-axis variable name and values
    x_axis_name = list(outputs.keys())[0] 
    x_axis_values = outputs[x_axis_name]
    # get y-axis variabe name and values --- output variable, usually Isp
    y_axis_name = list(outputs.keys())[1]
    y_axis_values = outputs[y_axis_name]
    # get apogee reached T/F and rocket parameters for each x,y pair
    apogee_reached_TF = outputs["apogee reached"]
    rocket_params = outputs["rocket parameters for given input"]
    
    # create main window
    root = tk.Tk()
    root.title("Parametric Analysis Dashboard")
    
    # create the matplotlib figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title(f'{graph_title}')
    ax.set_xlabel(f'{list(outputs.keys())[0]}')
    ax.set_ylabel(f'{list(outputs.keys())[1]}')
    
    # plot data
    data_point_colors = ["green" if reached else "red" for reached in apogee_reached_TF]
    ax.plot(x_axis_values, y_axis_values, color="black", linewidth=2, zorder=1)
    ax.grid(visible=True, which='major', color='lightgray', linestyle='-', linewidth=0.5)
    scatter = ax.scatter(x_axis_values, y_axis_values, c=data_point_colors, edgecolors="black", zorder=2)
    
    # Embed the matplotlib figure into the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    
    # Create a frame for rocket parameters
    param_frame = tk.Frame(root, bg="lightgray", padx=10, pady=10)
    param_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    param_label = tk.Label(param_frame, text="Rocket Parameters", bg="lightgray", font=("Arial", 14))
    param_label.pack(anchor="n", pady=10)

    param_text = tk.Text(param_frame, wrap=tk.WORD, bg="white", font=("Courier", 12), height=20, width=40)
    param_text.pack(fill=tk.BOTH, expand=True)

    # Function to display parameters when hovering over points
    param_text.insert("1.0", "Hover over data points to see rocket info.") # initial text to be displayed
    def on_hover(event):
        if event.inaxes == ax:
            for i, point in enumerate(scatter.get_offsets()):
                if event.xdata and event.ydata:
                    distance = ((point[0] - event.xdata)**2 + (point[1] - event.ydata)**2)**0.5
                    if distance < 0.2:  # adjust hover sensitivity
                        # clear and display the corresponding rocket parameters
                        param_text.delete("1.0", tk.END)
                        params = rocket_params[i]
                        for key, value in params.items():
                            param_text.insert(tk.END, f"{key}: {value}\n")
                        return

    # Connect the hover event to the canvas
    fig.canvas.mpl_connect("motion_notify_event", on_hover)

    # Start the Tkinter main loop
    root.mainloop()
    
    
# --------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY GRAPH: 2 INPUTS --------------------
# --------------------------------------------------------------------------

# plot a 3d graph of Isp vs two parametric inputs
def plot_two_inputs_vs_Isp(parametric_data):
    # Unpack data
    graph_title, inputs, outputs = unpack_data(parametric_data)
    # Extract x, y, and z data
    x_name = list(outputs.keys())[0]
    x_values = np.array(outputs[x_name])  # x-coordinates
    y_name = list(outputs.keys())[1]
    y_values = np.array(outputs[y_name])  # y-coordinates
    z_name = list(outputs.keys())[2]
    z_values = np.array(outputs[z_name])  # z-coordinates
    # Get apogee reached (True/False) for coloring
    apogee_reached_TF = outputs["apogee reached"]
    rocket_params = outputs["rocket parameters for given input"]
    
    # assign point colors
    data_point_colors = ["green" if reached else "red" for reached in apogee_reached_TF]
    
    # create figure with two panels for seeing the graph and the rocket parameters
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121, projection='3d')  # Left panel for the graph
    
    # Interpolate to create a smooth surface
    grid_x, grid_y = np.linspace(x_values.min(), x_values.max(), 50), np.linspace(y_values.min(), y_values.max(), 50)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)
    grid_z = griddata((x_values, y_values), z_values, (grid_x, grid_y), method='cubic')
    
    # Plot the surface with a red-to-green gradient along the z-axis
    truncated_brg = plt.cm.get_cmap("brg")(np.linspace(0.5, 1.0, 256)) # truncated section of brg so we get only the red-green gradient
    truncated_brg_cmap = ListedColormap(truncated_brg)
    surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=truncated_brg_cmap, alpha=0.6, edgecolor="none") # cmap is responsible for the color gradient. for more color gradient options see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    
    # Scatter points with colors based on apogee_reached
    scatter = ax.scatter(x_values, y_values, z_values, c=data_point_colors, s=50, edgecolors="black")
    
    # Add labels and title
    ax.set_xlabel(f"{x_name}")
    ax.set_ylabel(f"{y_name}")
    ax.set_zlabel(f"{z_name}")
    ax.set_title(f"{graph_title}")

    # Add a color bar
    #cbar = fig.colorbar(surf, ax=ax, pad=0.1)
    #cbar.set_label(z_name, rotation=270, labelpad=15)
    
    # Add text box on the right-hand side for dictionary details
    text_box_ax = fig.add_subplot(122)  # Right panel for the dictionary
    text_box_ax.axis("off")  # No axes for the text box
    text_box_ax.set_facecolor("lightgray") # doesn't do anything for some reason
    text_box = text_box_ax.text(0, 1, "Hover over a point to see details here.", verticalalignment="top", fontsize=10, family="monospace")
    
    # Add hover functionality to display dictionary content on the right
    cursor = mplcursors.cursor(scatter, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        idx = sel.index
        formatted_text = format_dict(rocket_params[idx])
        text_box.set_text(formatted_text)
        fig.canvas.draw_idle()

    # Function to format a dictionary neatly for display
    def format_dict(input_dictionary):
        lines = []
        max_key_length = max(len(str(key)) for key in input_dictionary.keys())
        for key, value in input_dictionary.items():
            dots = " " * (max_key_length - len(str(key)) + 0) #
            formatted_value = f"{value:.2f}" if isinstance(value, (float, int)) else value
            lines.append(f"{dots}{key}: {formatted_value}")
        return "\n".join(lines)

    # Display the popup window
    plt.show()
    
    
    
    
# --------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY GRAPH: 3 INPUTS --------------------
# --------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY GRAPH HELPER FUNCTIONS --------------------
# ---------------------------------------------------------------------------------


def unpack_data(parametric_data):
    graph_title = parametric_data["title"]
    inputs = parametric_data["parametric inputs"]
    outputs = parametric_data["parametric outputs"]
    return graph_title, inputs, outputs
