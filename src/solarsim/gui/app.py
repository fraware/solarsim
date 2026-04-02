"""Interactive 3D visualization and controls."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import (  # type: ignore[attr-defined]
    FigureCanvasTkAgg,
    NavigationToolbar2Tk,
)

from solarsim.body import Body
from solarsim.simulation import SolarSystemSim


def _traj_length(body: Body) -> int:
    return len(body.traj)


def _traj_array(body: Body) -> np.ndarray:
    return np.asarray(body.traj, dtype=np.float64)


class SolarSystemApp:
    """Tkinter shell around Matplotlib 3D and :class:`SolarSystemSim`."""

    def __init__(self, master: tk.Tk) -> None:
        self.master = master
        master.title("SolarSim — 3D solar system")

        try:
            scaling = master.winfo_fpixels("1i") / 72.0
            master.tk.call("tk", "scaling", scaling)
        except tk.TclError:
            pass

        self.sim = SolarSystemSim(base_dt=86400)
        self.paused = False

        self.dt_var = tk.DoubleVar(value=86400)
        self.speed_var = tk.IntVar(value=10)
        self.camera_mode = tk.StringVar(value="Free")
        self.follow_body = tk.StringVar(value="None")
        self.timeline_index = tk.IntVar(value=0)

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
        self.resume_button = tk.Button(ctrl_frame, text="Resume", command=self.resume_sim)
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

        self.time_label = tk.Label(ctrl_frame, text="Time: 0.0 days")
        self.time_label.pack(side=tk.LEFT, padx=10)
        self.energy_label = tk.Label(ctrl_frame, text="Energy: 0.0 J")
        self.energy_label.pack(side=tk.LEFT, padx=10)

        self.fig = plt.figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.ax.set_title("SolarSim — 3D")
        self.ax.set_xlabel("X (m)")
        self.ax.set_ylabel("Y (m)")
        self.ax.set_zlabel("Z (m)")
        self.ax.set_facecolor("black")
        self.fig.patch.set_facecolor("black")
        self.set_default_axes()

        np.random.seed(42)
        stars = np.random.uniform(-6e12, 6e12, (3, 300))
        self.ax.scatter(stars[0], stars[1], stars[2], color="white", s=1, alpha=0.3)

        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, master)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)

        self.body_plots: dict[str, Any] = {}
        self.traj_plots: dict[str, Any] = {}
        self.body_labels: dict[str, Any] = {}
        for body in self.sim.bodies:
            (point,) = self.ax.plot([], [], [], "o", color=body.color, label=body.name)
            (line,) = self.ax.plot([], [], [], "-", color=body.color, linewidth=0.5)
            self.body_plots[body.name] = point
            self.traj_plots[body.name] = line
            lbl = self.ax.text(
                float(body.pos[0]),
                float(body.pos[1]),
                float(body.pos[2]),
                body.name,
                color=body.color,
                fontsize=8,
            )
            self.body_labels[body.name] = lbl
        self.ax.legend(loc="upper right", fontsize="small")

        self.ani = FuncAnimation(self.fig, self.update, interval=50, blit=False)

    def set_default_axes(self) -> None:
        self.ax.set_xlim(-6e12, 6e12)
        self.ax.set_ylim(-6e12, 6e12)
        self.ax.set_zlim(-6e12, 6e12)

    def _max_frame_index(self) -> int:
        if not self.sim.bodies:
            return 0
        return max(_traj_length(b) for b in self.sim.bodies) - 1

    def _remove_plot_artists_for_destroyed_bodies(self) -> None:
        alive = {b.name for b in self.sim.bodies}
        for name in list(self.body_plots):
            if name not in alive:
                self.body_plots[name].remove()
                self.traj_plots[name].remove()
                self.body_labels[name].remove()
                del self.body_plots[name]
                del self.traj_plots[name]
                del self.body_labels[name]

    def update(self, frame: int) -> list[Any]:
        if not self.paused:
            steps = self.speed_var.get()
            self.sim.base_dt = self.dt_var.get()
            for _ in range(steps):
                self.sim.step()
            current_frame = self._max_frame_index()
            self.timeline_slider.config(to=max(current_frame, 0))
            self.timeline_index.set(current_frame)

        self._remove_plot_artists_for_destroyed_bodies()

        idx = self.timeline_index.get()
        for body in self.sim.bodies:
            traj = _traj_array(body)
            pos = traj[idx] if idx < len(traj) else body.pos
            self.body_plots[body.name].set_data((pos[0],), (pos[1],))
            self.body_plots[body.name].set_3d_properties(pos[2])
            end = min(idx + 1, len(traj))
            self.traj_plots[body.name].set_data(traj[:end, 0], traj[:end, 1])
            self.traj_plots[body.name].set_3d_properties(traj[:end, 2])
            lbl = self.body_labels[body.name]
            lbl.set_position((float(pos[0]), float(pos[1])))
            lbl.set_3d_properties(float(pos[2]))

        days = self.sim.time / 86400.0
        self.time_label.config(text=f"Time: {days:.2f} days")
        energy = self.sim.compute_total_energy()
        self.energy_label.config(text=f"Energy: {energy:.2e} J")

        mode = self.camera_mode.get()
        if mode == "Follow" and self.follow_body.get() != "None":
            body_name = self.follow_body.get()
            b = next((x for x in self.sim.bodies if x.name == body_name), None)
            if b is not None:
                traj = _traj_array(b)
                cpos = traj[idx] if idx < len(traj) else b.pos
                delta = 1e11
                self.ax.set_xlim(cpos[0] - delta, cpos[0] + delta)
                self.ax.set_ylim(cpos[1] - delta, cpos[1] + delta)
                self.ax.set_zlim(cpos[2] - delta, cpos[2] + delta)
        elif mode == "Top-Down":
            self.ax.view_init(elev=90, azim=-90)
        else:
            self.set_default_axes()

        self.canvas.draw()
        return []

    def pause_sim(self) -> None:
        self.paused = True

    def resume_sim(self) -> None:
        self.paused = False

    def step_sim(self) -> None:
        if self.paused:
            self.sim.step()
            current_frame = self._max_frame_index()
            self.timeline_slider.config(to=max(current_frame, 0))
            self.timeline_index.set(current_frame)
            _ = self.update(0)


def run_app() -> None:
    root = tk.Tk()
    SolarSystemApp(root)
    root.mainloop()
