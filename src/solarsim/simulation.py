"""N-body integration: velocity Verlet, adaptive timestep, collisions, energy."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

from solarsim.body import Body
from solarsim.constants import C_LIGHT, G
from solarsim.ephemeris import create_default_bodies
from solarsim.physics_numba import accelerations_kernel, stack_state


class SolarSystemSim:
    """N-body solar model: velocity Verlet, adaptive timestep, optional collisions.

    Gravity uses Plummer softening ``sqrt(r^2 + eps^2)`` to limit singularity spikes.
    For the body named ``Mercury``, an extra factor ``1 + 3 l^2 / (r^2 c^2)`` is applied
    to the Newtonian acceleration (``l`` = specific angular momentum). This mimics the
    *order of magnitude* of GR perihelion effects for teaching; it is not a valid
    substitute for relativistic orbit propagation.
    """

    def __init__(
        self,
        base_dt: float,
        *,
        max_traj_points: int | None = 100_000,
        softening_m: float = 1e6,
        body_factory: Callable[[], list[Body]] | None = None,
    ) -> None:
        self.base_dt = base_dt
        self.dt = base_dt
        self.time = 0.0
        self._max_traj_points = max_traj_points
        self._eps2 = float(softening_m) ** 2
        factory = body_factory or (lambda: create_default_bodies(max_traj_points=max_traj_points))
        self.bodies: list[Body] = factory()
        self._mercury_idx = next(
            (i for i, b in enumerate(self.bodies) if b.name == "Mercury"),
            -1,
        )

    def compute_accelerations(self) -> list[np.ndarray]:
        pos, vel, mass = stack_state(self.bodies)
        acc = accelerations_kernel(pos, vel, mass, self._mercury_idx, G, C_LIGHT, self._eps2)
        return [acc[i].copy() for i in range(len(self.bodies))]

    def velocity_verlet_step(self) -> None:
        accs = self.compute_accelerations()
        dt = self.dt
        for i, body in enumerate(self.bodies):
            body.pos = body.pos + body.vel * dt + 0.5 * accs[i] * dt**2
        new_accs = self.compute_accelerations()
        for i, body in enumerate(self.bodies):
            body.vel = body.vel + 0.5 * (accs[i] + new_accs[i]) * dt
            body.append_traj(body.pos)
        self.time += dt

    def adaptive_dt(self) -> None:
        min_dist = math.inf
        n = len(self.bodies)
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(self.bodies[j].pos - self.bodies[i].pos))
                if d < min_dist:
                    min_dist = d
        d_min = 1e9
        d_max = 1e11
        factor = 1.0
        if min_dist < d_min:
            factor = 0.5
        elif min_dist > d_max:
            factor = 1.1
        self.dt = float(np.clip(self.dt * factor, self.base_dt / 10, self.base_dt * 10))

    def detect_collisions(self) -> None:
        i = 0
        while i < len(self.bodies):
            j = i + 1
            while j < len(self.bodies):
                d = float(np.linalg.norm(self.bodies[i].pos - self.bodies[j].pos))
                if d < (self.bodies[i].radius + self.bodies[j].radius):
                    if self.bodies[i].mass >= self.bodies[j].mass:
                        larger, smaller = self.bodies[i], self.bodies[j]
                        remove_idx = j
                    else:
                        larger, smaller = self.bodies[j], self.bodies[i]
                        remove_idx = i
                    total_mass = larger.mass + smaller.mass
                    new_vel = (larger.mass * larger.vel + smaller.mass * smaller.vel) / total_mass
                    new_pos = (larger.mass * larger.pos + smaller.mass * smaller.pos) / total_mass
                    larger.mass = total_mass
                    larger.vel = new_vel
                    larger.pos = new_pos
                    larger.radius = (larger.radius**3 + smaller.radius**3) ** (1 / 3)
                    larger.extend_traj_from(smaller)
                    del self.bodies[remove_idx]
                    if remove_idx == i:
                        j = i + 1
                    continue
                j += 1
            i += 1
        self._mercury_idx = next(
            (k for k, b in enumerate(self.bodies) if b.name == "Mercury"),
            -1,
        )

    def compute_total_energy(self) -> float:
        kinetic = 0.0
        potential = 0.0
        n = len(self.bodies)
        for body in self.bodies:
            kinetic += 0.5 * body.mass * float(np.dot(body.vel, body.vel))
        for i in range(n):
            for j in range(i + 1, n):
                r_vec = self.bodies[j].pos - self.bodies[i].pos
                r = math.sqrt(float(np.dot(r_vec, r_vec)) + self._eps2)
                potential -= G * self.bodies[i].mass * self.bodies[j].mass / r
        return kinetic + potential

    def step(self) -> None:
        self.velocity_verlet_step()
        self.adaptive_dt()
        self.detect_collisions()
