import matplotlib.pyplot as plt

def graph_everything(timesteps_dict):
    """Plot position, velocity, acceleration, and forces in separate windows simultaneously."""
    plot_position(timesteps_dict, block=False)
    plot_velocity(timesteps_dict, block=False)
    plot_acceleration(timesteps_dict, block=False)
    plot_forces(timesteps_dict, block=False)
    plt.show()  # Keeps all windows open until manually closed

def plot_position(timesteps_dict, block=True):
    """Plot position as a function of time."""
    time = timesteps_dict["time"]
    position = timesteps_dict["position"]
    plt.figure()
    plt.plot(time, position, label="Position", linewidth=2)
    plt.title("Rocket Position vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Position (m)")
    plt.grid(True)
    plt.legend()
    plt.show(block=block)

def plot_velocity(timesteps_dict, block=True):
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
    plt.show(block=block)

def plot_acceleration(timesteps_dict, block=True):
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
    plt.show(block=block)

def plot_forces(timesteps_dict, block=True):
    """Plot net force, drag force, and gravitational force as functions of time."""
    time = timesteps_dict["time"]
    net_force = timesteps_dict["net force"]
    drag_force = timesteps_dict["drag force"]
    gravitational_force = timesteps_dict["gravitational force"]

    plt.figure()
    plt.plot(time, net_force, label="Net Force", linewidth=2, color="blue")
    plt.plot(time, drag_force, label="Drag Force", linewidth=2, color="red")
    plt.plot(time, gravitational_force, label="Gravitational Force", linewidth=2, color="purple")
    plt.title("Rocket Forces vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.grid(True)
    plt.legend()
    plt.show(block=block)
