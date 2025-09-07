import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
import numpy as np
from matplotlib.widgets import Slider
from matplotlib.patches import Rectangle




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


# -------------------------------------------------------------------
# -------------------- PARAMETRIC STUDY PLOTTING --------------------
# -------------------------------------------------------------------

def plot_one_input_vs_Isp(graph_dict):
    """Plot the graph with dynamic x and y labels, color-coded data points, connected lines, and hoverable information."""
    # Extract keys and values from the dictionary
    x_key = [key for key in graph_dict.keys() if "oxidizer" in key.lower()][0]  # Dynamic x-axis key
    y_key = "specific impulse"  # Static y-axis key
    x_label = x_key.capitalize()
    y_label = y_key.capitalize()

    # Extract data
    x_values = graph_dict[x_key]
    y_values = graph_dict[y_key]
    apogee_reached = graph_dict["apogee reached"]
    rocket_params = graph_dict["rocket parameters for given input"]

    # Create a figure with two panels
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(1, 2, width_ratios=[3, 1])  # 2 panels, graph takes up more space
    ax_plot = fig.add_subplot(gs[0])
    ax_dict = fig.add_subplot(gs[1])

    # Plot the graph
    colors = ['green' if reached else 'red' for reached in apogee_reached]
    ax_plot.plot(x_values, y_values, linestyle='-', color='blue', alpha=0.5)  # Connect the points with a line
    scatter = ax_plot.scatter(x_values, y_values, c=colors, label=f'{y_label} vs {x_label}', zorder=3)

    # Add labels, title, and legend
    ax_plot.set_xlabel(x_label)
    ax_plot.set_ylabel(y_label)
    ax_plot.set_title(f'{y_label} vs {x_label}')
    ax_plot.legend()
    ax_plot.grid(True)

    # Initialize the dictionary display
    create_scrollable_container(ax_dict, {"Hover over a point to see details": ""})
    enable_scrolling(ax_dict, fig, rocket_params[0])

    # Add hover functionality
    annot = ax_plot.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        """Update annotation with the corresponding rocket parameter dictionary."""
        idx = ind["ind"][0]
        pos = scatter.get_offsets()[idx]

        # Check if the popup is within bounds, adjust if necessary
        x_pos, y_pos = pos
        bbox = annot.get_window_extent(renderer=fig.canvas.get_renderer())
        if bbox.x1 > ax_plot.get_window_extent().x1:  # Overflow to the right
            annot.xy = (x_pos - 0.1, y_pos - 0.1)
        elif bbox.y0 < ax_plot.get_window_extent().y0:  # Overflow at the bottom
            annot.xy = (x_pos - 0.1, y_pos + 0.1)
        else:
            annot.xy = pos

        text = f"{x_label}: {x_values[idx]:.2f}\n{y_label}: {y_values[idx]:.2f}"
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(colors[idx])
        annot.get_bbox_patch().set_alpha(0.8)

        # Update the right-hand panel with the dictionary
        create_scrollable_container(ax_dict, rocket_params[idx])
        fig.canvas.draw_idle()

    def hover(event):
        """Handle hover events on the graph."""
        vis = annot.get_visible()
        if event.inaxes == ax_plot:
            cont, ind = scatter.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.tight_layout()
    plt.show()







# takes 3 lists (each list must have the same amount of items), outputs a 3d graph
def plot_OF_vs_CP_vs_Isp(of_ratio, chamber_pressure, specific_impulse):    
    # create 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # make it look nice
    #ax.scatter(of_ratio, chamber_pressure, specific_impulse, c=specific_impulse, cmap='viridis', marker='o')
    ax.plot_trisurf(of_ratio, chamber_pressure, specific_impulse, cmap='viridis', edgecolor='none')
    
    # set axis labels
    ax.set_xlabel('OF Ratio')
    ax.set_ylabel('Chamber Pressure (psi)')
    ax.set_zlabel('Specific Impulse (s)')
    
    ax.set_title('Specific Impulse as a function of OF ratio and chamber pressure')
    
def plot_Mo_vs_Lf_vs_Isp(ox_dot, fuel_length, specific_impulse):
    # create 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot_trisurf(ox_dot, fuel_length, specific_impulse, cmap='viridis', edgecolor='none')
    
    # set axis labels
    ax.set_xlabel('Oxidizer mass flow rate (kg/s)')
    ax.set_ylabel('Fuel length (in)')
    ax.set_zlabel('Specific Impulse (s)')
    
    ax.set_title('Specific Impulse as a function of ox flow rate and fuel length')
    
# ideal rocket nozzle radius as a function of oxidizer mass flow rate and fuel length
def plot_Mo_vs_Lf_vs_Re(ox_dot, fuel_length, nozzle_outlet_radius):
    # create 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot_trisurf(ox_dot, fuel_length, nozzle_outlet_radius, cmap='viridis', edgecolor='none')
    
    # set axis labels
    ax.set_xlabel('Oxidizer mass flow rate (kg/s)')
    ax.set_ylabel('Fuel length (in)')
    ax.set_zlabel('Nozzle outlet radius (in)')
    
    ax.set_title('Nozzle outlet radius as a function of ox flow rate and fuel length')
    
    
    
    
    
    
    
    
    
    
    
# helper function to print dictionaries nicely
def print_dict_in_panel(ax, input_dictionary):
    """Helper function to display a dictionary in the provided subplot with scrolling."""
    ax.clear()  # Clear the panel
    ax.axis("off")  # Turn off the axis
    y_pos = 1.0  # Start at the top of the panel
    line_spacing = 0.05  # Space between lines

    # Format the dictionary values
    formatted_dict = {key: format_value(value) for key, value in input_dictionary.items()}

    # Display each key-value pair with scrolling enabled
    ax._scrollable_text = [
        f"{key}: {value}" for key, value in formatted_dict.items()
    ]
    max_visible_lines = 20
    current_view_offset = getattr(ax, "_scroll_offset", 0)

    # Draw the visible lines
    for i, line in enumerate(ax._scrollable_text[current_view_offset:current_view_offset + max_visible_lines]):
        ax.text(0, y_pos, line, fontsize=10, va="top", ha="left")
        y_pos -= line_spacing

def format_value(value):
    """Round numerical values to 2 decimals or return as-is if non-numeric."""
    if isinstance(value, (int, float, np.floating)):
        return round(value, 2)
    return value

def enable_scrolling(ax, fig):
    """Add scrolling functionality to the dictionary display panel."""
    ax._scroll_offset = 0  # Initialize scroll offset

    def on_scroll(event):
        if event.inaxes == ax:
            scroll_step = 1  # Scroll by one line
            if event.button == "up":
                ax._scroll_offset = max(ax._scroll_offset - scroll_step, 0)
            elif event.button == "down":
                max_offset = max(len(ax._scrollable_text) - 20, 0)
                ax._scroll_offset = min(ax._scroll_offset + scroll_step, max_offset)

            print_dict_in_panel(ax, ax._last_dictionary)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("scroll_event", on_scroll)
    
def format_value(value):
    """Round numerical values to 2 decimals or return as-is if non-numeric."""
    if isinstance(value, (int, float, np.floating)):
        return round(value, 2)
    return value


def create_scrollable_container(ax, input_dictionary):
    """Helper function to create a scrollable container for the dictionary."""
    ax.clear()  # Clear the panel
    ax.axis("off")  # Turn off the axis
    
    # Add a light gray background
    ax.add_patch(Rectangle((0, 0), 1, 1, color='lightgray', transform=ax.transAxes, zorder=0))
    
    # Add a bolded header
    ax.text(0.02, 1.0, "Rocket Parameters", fontsize=12, weight='bold', ha='left', va='top', transform=ax.transAxes)
    
    # Format the dictionary values
    formatted_dict = {key: format_value(value) for key, value in input_dictionary.items()}
    
    # Scrollable content
    scroll_text = [
        f"{key}: {value}" for key, value in formatted_dict.items()
    ]
    max_visible_lines = 20
    scroll_offset = getattr(ax, "_scroll_offset", 0)
    
    # Display content within scrollable container
    y_pos = 0.95  # Start just below the header
    line_spacing = 0.035  # Space between lines
    for i, line in enumerate(scroll_text[scroll_offset:scroll_offset + max_visible_lines]):
        ax.text(0.02, y_pos, line, fontsize=10, ha="left", va="top", transform=ax.transAxes)
        y_pos -= line_spacing


def enable_scrolling(ax, fig, input_dictionary):
    """Add scrolling functionality to the dictionary display panel."""
    ax._scroll_offset = 0  # Initialize scroll offset

    def on_scroll(event):
        if event.inaxes == ax:
            scroll_step = 1  # Scroll by one line
            if event.button == "up":
                ax._scroll_offset = max(ax._scroll_offset - scroll_step, 0)
            elif event.button == "down":
                max_offset = max(len(input_dictionary) - 20, 0)
                ax._scroll_offset = min(ax._scroll_offset + scroll_step, max_offset)

            create_scrollable_container(ax, input_dictionary)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("scroll_event", on_scroll)