#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.mplot3d import Axes3D  # enable 3D plotting
from numba import njit

# Physical constants
G = 6.67430e-11  # gravitational constant [m^3 kg^-1 s^-2]
c = 3e8  # speed of light [m/s]


# =============================================================================
# Body class – represents a celestial body.
# =============================================================================
class Body:
    def __init__(self, name, mass, pos, vel, color, radius):
        self.name = name
        self.mass = mass
        self.pos = np.array(pos, dtype=float)  # 3D position vector
        self.vel = np.array(vel, dtype=float)  # 3D velocity vector
        self.color = color
        self.radius = radius  # collision radius [m]
        self.traj = [self.pos.copy()]  # trajectory history


# =============================================================================
# SolarSystemSim – holds the simulation state and performs integration.
# =============================================================================
class SolarSystemSim:
    def __init__(self, base_dt):
        self.base_dt = base_dt  # base time step from UI (seconds)
        self.dt = base_dt  # current adaptive time step
        self.time = 0.0  # simulation time in seconds
        self.bodies = []
        self.init_bodies()

    def init_bodies(self):
        # Here we use approximate data. In a real version, you might import
        # precise orbital elements from AstroPy.
        # For collision radii, we use rough approximations.
        # (Units: mass in kg, distances in m, velocity in m/s, radius in m)
        self.bodies.append(Body("Sun", 1.989e30, [0, 0, 0], [0, 0, 0], "yellow", 7e8))
        self.bodies.append(
            Body(
                "Mercury",
                3.3e23,
                [5.79e10, 0, 0],
                [0, 47400 * np.cos(np.deg2rad(7)), 47400 * np.sin(np.deg2rad(7))],
                "gray",
                2.4e6,
            )
        )
        self.bodies.append(
            Body(
                "Venus",
                4.87e24,
                [1.082e11, 0, 0],
                [0, 35000 * np.cos(np.deg2rad(3.4)), 35000 * np.sin(np.deg2rad(3.4))],
                "orange",
                6e6,
            )
        )
        self.bodies.append(
            Body("Earth", 5.97e24, [1.496e11, 0, 0], [0, 29780, 0], "blue", 6.4e6)
        )
        # Earth’s Moon (placed relative to Earth)
        moon_dist = 3.84e8
        self.bodies.append(
            Body(
                "Moon",
                7.35e22,
                [1.496e11 + moon_dist, 0, 0],
                [0, 29780 + 1022, 0],
                "lightgray",
                1.7e6,
            )
        )
        self.bodies.append(
            Body(
                "Mars",
                6.39e23,
                [2.279e11, 0, 0],
                [0, 24077 * np.cos(np.deg2rad(1.85)), 24077 * np.sin(np.deg2rad(1.85))],
                "red",
                3.4e6,
            )
        )
        self.bodies.append(
            Body(
                "Jupiter",
                1.898e27,
                [7.78e11, 0, 0],
                [0, 13070 * np.cos(np.deg2rad(1.3)), 13070 * np.sin(np.deg2rad(1.3))],
                "saddlebrown",
                7e7,
            )
        )
        self.bodies.append(
            Body(
                "Saturn",
                5.683e26,
                [1.43e12, 0, 0],
                [0, 9690 * np.cos(np.deg2rad(2.5)), 9690 * np.sin(np.deg2rad(2.5))],
                "goldenrod",
                6e7,
            )
        )
        self.bodies.append(
            Body(
                "Uranus",
                8.681e25,
                [2.87e12, 0, 0],
                [0, 6810 * np.cos(np.deg2rad(0.77)), 6810 * np.sin(np.deg2rad(0.77))],
                "lightblue",
                2.5e7,
            )
        )
        self.bodies.append(
            Body(
                "Neptune",
                1.024e26,
                [4.5e12, 0, 0],
                [0, 5430 * np.cos(np.deg2rad(1.77)), 5430 * np.sin(np.deg2rad(1.77))],
                "darkblue",
                2.5e7,
            )
        )
        self.bodies.append(
            Body(
                "Pluto",
                1.309e22,
                [5.9e12, 0, 0],
                [0, 4740 * np.cos(np.deg2rad(17)), 4740 * np.sin(np.deg2rad(17))],
                "purple",
                1.2e6,
            )
        )

    def compute_accelerations(self):
        # Compute Newtonian accelerations for each body.
        n = len(self.bodies)
        accs = [np.zeros(3) for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                r_vec = self.bodies[j].pos - self.bodies[i].pos
                r = np.linalg.norm(r_vec)
                if r == 0:
                    continue
                # Newtonian acceleration magnitude:
                a_newt = G * self.bodies[j].mass / r**3
                # Basic acceleration contribution:
                a = a_newt * r_vec
                # If body i is Mercury, add a simple relativistic correction:
                if self.bodies[i].name == "Mercury":
                    # Compute specific angular momentum
                    l_vec = np.cross(self.bodies[i].pos, self.bodies[i].vel)
                    l2 = np.dot(l_vec, l_vec)
                    rel_corr = 1 + 3 * l2 / (r**2 * c**2)
                    a *= rel_corr
                accs[i] += a
        return accs

    def velocity_verlet_step(self):
        # Perform one integration step using velocity Verlet
        accs = self.compute_accelerations()
        # First, update positions:
        for i, body in enumerate(self.bodies):
            body.pos = body.pos + body.vel * self.dt + 0.5 * accs[i] * self.dt**2
        # Compute new accelerations at updated positions:
        new_accs = self.compute_accelerations()
        # Update velocities:
        for i, body in enumerate(self.bodies):
            body.vel = body.vel + 0.5 * (accs[i] + new_accs[i]) * self.dt
            body.traj.append(body.pos.copy())
        self.time += self.dt

    def adaptive_dt(self):
        # Adjust dt based on minimum distance between any two bodies.
        min_dist = np.inf
        n = len(self.bodies)
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(self.bodies[j].pos - self.bodies[i].pos)
                if d < min_dist:
                    min_dist = d
        # Define thresholds (these numbers can be tuned)
        d_min = 1e9  # if too close, reduce dt
        d_max = 1e11  # if very far, we can increase dt (up to base_dt)
        # Simple adaptive scheme:
        factor = 1.0
        if min_dist < d_min:
            factor = 0.5
        elif min_dist > d_max:
            factor = 1.1
        # Clamp dt between 1/10th and 10x the base dt:
        new_dt = np.clip(self.dt * factor, self.base_dt / 10, self.base_dt * 10)
        self.dt = new_dt

    def detect_collisions(self):
        # Check for collisions between bodies.
        # If two bodies come within the sum of their radii, merge them.
        i = 0
        while i < len(self.bodies):
            j = i + 1
            while j < len(self.bodies):
                d = np.linalg.norm(self.bodies[i].pos - self.bodies[j].pos)
                if d < (self.bodies[i].radius + self.bodies[j].radius):
                    # Merge smaller into larger.
                    if self.bodies[i].mass >= self.bodies[j].mass:
                        larger, smaller = self.bodies[i], self.bodies[j]
                    else:
                        larger, smaller = self.bodies[j], self.bodies[i]
                    # Conservation of mass and momentum:
                    total_mass = larger.mass + smaller.mass
                    new_vel = (
                        larger.mass * larger.vel + smaller.mass * smaller.vel
                    ) / total_mass
                    # Weighted position:
                    new_pos = (
                        larger.mass * larger.pos + smaller.mass * smaller.pos
                    ) / total_mass
                    larger.mass = total_mass
                    larger.vel = new_vel
                    larger.pos = new_pos
                    larger.radius = (larger.radius**3 + smaller.radius**3) ** (
                        1 / 3
                    )  # approximate new radius
                    # Append smaller's trajectory to larger's (for visualization)
                    larger.traj += smaller.traj
                    # Remove the smaller body:
                    del self.bodies[j]
                else:
                    j += 1
            i += 1

    def compute_total_energy(self):
        # Compute total energy (kinetic + potential) for monitoring conservation.
        kinetic = 0.0
        potential = 0.0
        n = len(self.bodies)
        for body in self.bodies:
            kinetic += 0.5 * body.mass * np.dot(body.vel, body.vel)
        for i in range(n):
            for j in range(i + 1, n):
                r = np.linalg.norm(self.bodies[j].pos - self.bodies[i].pos)
                potential -= G * self.bodies[i].mass * self.bodies[j].mass / r
        return kinetic + potential

    def step(self):
        # One simulation step: perform integration, adaptive dt, and collision detection.
        self.velocity_verlet_step()
        self.adaptive_dt()
        self.detect_collisions()


# =============================================================================
# Tkinter application with interactive controls and enhanced 3D visualization.
# =============================================================================
class SolarSystemApp:
    def __init__(self, master):
        self.master = master
        master.title("Enhanced Solar System 3D Simulation")

        # Simulation instance with base dt (default 1 day)
        self.sim = SolarSystemSim(base_dt=86400)
        self.paused = False
        self.use_gpu = tk.BooleanVar(value=False)  # Placeholder for GPU mode

        # Control variables
        self.dt_var = tk.DoubleVar(value=86400)
        self.speed_var = tk.IntVar(value=10)  # steps per frame
        self.camera_mode = tk.StringVar(value="Free")
        self.follow_body = tk.StringVar(value="None")
        self.timeline_index = tk.IntVar(value=0)

        # Setup control panel
        ctrl_frame = tk.Frame(master)
        ctrl_frame.pack(side=tk.TOP, fill=tk.X)

        tk.Label(ctrl_frame, text="Base Time Step (s):").pack(side=tk.LEFT)
        self.dt_entry = tk.Entry(ctrl_frame, textvariable=self.dt_var, width=10)
        self.dt_entry.pack(side=tk.LEFT)

        tk.Label(ctrl_frame, text="Steps/Frame:").pack(side=tk.LEFT)
        self.speed_slider = tk.Scale(
            ctrl_frame, variable=self.speed_var, from_=1, to=100, orient=tk.HORIZONTAL
        )
        self.speed_slider.pack(side=tk.LEFT)

        self.pause_button = tk.Button(ctrl_frame, text="Pause", command=self.pause_sim)
        self.pause_button.pack(side=tk.LEFT, padx=5)
        self.resume_button = tk.Button(
            ctrl_frame, text="Resume", command=self.resume_sim
        )
        self.resume_button.pack(side=tk.LEFT, padx=5)
        self.step_button = tk.Button(ctrl_frame, text="Step", command=self.step_sim)
        self.step_button.pack(side=tk.LEFT, padx=5)

        tk.Label(ctrl_frame, text="Camera:").pack(side=tk.LEFT, padx=5)
        cam_options = ["Free", "Follow", "Top-Down"]
        self.cam_menu = ttk.Combobox(
            ctrl_frame,
            textvariable=self.camera_mode,
            values=cam_options,
            state="readonly",
            width=10,
        )
        self.cam_menu.pack(side=tk.LEFT)

        tk.Label(ctrl_frame, text="Follow:").pack(side=tk.LEFT, padx=5)
        follow_options = ["None"] + [body.name for body in self.sim.bodies]
        self.follow_menu = ttk.Combobox(
            ctrl_frame,
            textvariable=self.follow_body,
            values=follow_options,
            state="readonly",
            width=10,
        )
        self.follow_menu.pack(side=tk.LEFT)

        # GPU acceleration toggle (placeholder)
        self.gpu_chk = tk.Checkbutton(
            ctrl_frame, text="GPU Acceleration", variable=self.use_gpu
        )
        self.gpu_chk.pack(side=tk.LEFT, padx=5)

        # Timeline slider for rewinding/scrubbing
        self.timeline_slider = tk.Scale(
            ctrl_frame,
            variable=self.timeline_index,
            from_=0,
            to=0,
            orient=tk.HORIZONTAL,
            label="Timeline",
            length=200,
        )
        self.timeline_slider.pack(side=tk.LEFT, padx=10)

        # Data overlay labels
        self.time_label = tk.Label(ctrl_frame, text="Time: 0.0 days")
        self.time_label.pack(side=tk.LEFT, padx=10)
        self.energy_label = tk.Label(ctrl_frame, text="Energy: 0.0 J")
        self.energy_label.pack(side=tk.LEFT, padx=10)

        # Setup Matplotlib 3D figure
        self.fig = plt.figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.ax.set_title("Enhanced Solar System 3D Simulation")
        self.ax.set_xlabel("X (m)")
        self.ax.set_ylabel("Y (m)")
        self.ax.set_zlabel("Z (m)")
        self.ax.set_facecolor("black")
        self.fig.patch.set_facecolor("black")
        self.set_default_axes()

        # Add background stars
        np.random.seed(42)
        stars = np.random.uniform(-6e12, 6e12, (3, 300))
        self.ax.scatter(stars[0], stars[1], stars[2], color="white", s=1, alpha=0.3)

        # Embed Matplotlib canvas into Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, master)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)

        # Create plot objects for bodies and trajectories
        self.body_plots = {}
        self.traj_plots = {}
        for body in self.sim.bodies:
            (point,) = self.ax.plot([], [], [], "o", color=body.color, label=body.name)
            (line,) = self.ax.plot([], [], [], "-", color=body.color, linewidth=0.5)
            self.body_plots[body.name] = point
            self.traj_plots[body.name] = line
            # Also add a label near each body
            self.ax.text(
                body.pos[0], body.pos[1], body.pos[2], body.name, color=body.color
            )
        self.ax.legend(loc="upper right", fontsize="small")

        # Animation: update every 50 ms
        self.ani = FuncAnimation(self.fig, self.update, interval=50, blit=False)

    def set_default_axes(self):
        self.ax.set_xlim(-6e12, 6e12)
        self.ax.set_ylim(-6e12, 6e12)
        self.ax.set_zlim(-6e12, 6e12)

    def update(self, frame):
        # If not scrubbing via timeline, advance simulation:
        if not self.paused:
            steps = self.speed_var.get()
            # Update base_dt from control panel
            self.sim.base_dt = self.dt_var.get()
            for _ in range(steps):
                self.sim.step()
            # Update timeline slider maximum
            current_frame = len(self.sim.bodies[0].traj) - 1
            self.timeline_slider.config(to=current_frame)
            self.timeline_index.set(current_frame)

        # Determine which frame to display (timeline slider value)
        idx = self.timeline_index.get()
        # Update each body's plot from its trajectory history:
        for body in self.sim.bodies:
            traj = np.array(body.traj)
            if idx < len(traj):
                pos = traj[idx]
            else:
                pos = body.pos
            self.body_plots[body.name].set_data(pos[0], pos[1])
            self.body_plots[body.name].set_3d_properties(pos[2])
            # Update trajectory line (show entire history up to idx)
            self.traj_plots[body.name].set_data(traj[:idx, 0], traj[:idx, 1])
            self.traj_plots[body.name].set_3d_properties(traj[:idx, 2])

        # Data overlay: simulation time in days and total energy
        days = self.sim.time / 86400.0
        self.time_label.config(text=f"Time: {days:.2f} days")
        energy = self.sim.compute_total_energy()
        self.energy_label.config(text=f"Energy: {energy:.2e} J")

        # Camera control
        mode = self.camera_mode.get()
        if mode == "Follow" and self.follow_body.get() != "None":
            # Center view on the selected body
            body_name = self.follow_body.get()
            b = next((body for body in self.sim.bodies if body.name == body_name), None)
            if b is not None:
                delta = 1e11
                self.ax.set_xlim(b.pos[0] - delta, b.pos[0] + delta)
                self.ax.set_ylim(b.pos[1] - delta, b.pos[1] + delta)
                self.ax.set_zlim(b.pos[2] - delta, b.pos[2] + delta)
        elif mode == "Top-Down":
            # Top-down view: fix the view angles
            self.ax.view_init(elev=90, azim=-90)
        else:
            # Free view: reset to default limits if not following a body
            self.set_default_axes()

        self.canvas.draw()

    def pause_sim(self):
        self.paused = True

    def resume_sim(self):
        self.paused = False

    def step_sim(self):
        if self.paused:
            self.sim.step()
            current_frame = len(self.sim.bodies[0].traj) - 1
            self.timeline_slider.config(to=current_frame)
            self.timeline_index.set(current_frame)
            self.update(0)


if __name__ == "__main__":
    root = tk.Tk()
    app = SolarSystemApp(root)
    root.mainloop()
