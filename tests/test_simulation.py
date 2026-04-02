"""Headless tests for integration, collisions, and timestep policy."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import pytest

from solarsim.body import Body
from solarsim.constants import G
from solarsim.simulation import SolarSystemSim


def _two_approaching() -> list[Body]:
    return [
        Body("A", 1.0, [-10.0, 0, 0], [5.0, 0, 0], "red", 1.0, max_traj_points=None),
        Body("B", 1.0, [10.0, 0, 0], [-5.0, 0, 0], "blue", 1.0, max_traj_points=None),
    ]


def test_collision_conserves_momentum_and_mass() -> None:
    factory: Callable[[], list[Body]] = _two_approaching
    sim = SolarSystemSim(
        0.05,
        max_traj_points=None,
        softening_m=1e-6,
        body_factory=factory,
    )
    for _ in range(500):
        sim.step()
        if len(sim.bodies) == 1:
            break
    assert len(sim.bodies) == 1
    b = sim.bodies[0]
    assert b.mass == pytest.approx(2.0)
    assert b.vel[0] == pytest.approx(0.0, abs=1e-6)


def test_adaptive_dt_stays_within_clamp() -> None:
    sim = SolarSystemSim(
        100.0,
        max_traj_points=None,
        softening_m=1e6,
        body_factory=_two_approaching,
    )
    for _ in range(50):
        sim.step()
        assert sim.dt >= sim.base_dt / 10 - 1e-9
        assert sim.dt <= sim.base_dt * 10 + 1e-9


def test_circular_orbit_radius_stable() -> None:
    def factory() -> list[Body]:
        m_central = 1.989e30
        r = 1.496e11
        v_circ = math.sqrt(G * m_central / r)
        return [
            Body("Sun", m_central, [0, 0, 0], [0, 0, 0], "yellow", 7e8, max_traj_points=5000),
            Body(
                "Probe",
                1.0,
                [r, 0, 0],
                [0, v_circ, 0],
                "cyan",
                1e3,
                max_traj_points=5000,
            ),
        ]

    sim = SolarSystemSim(
        3600.0,
        max_traj_points=5000,
        softening_m=1e8,
        body_factory=factory,
    )
    r0 = float(np.linalg.norm(sim.bodies[1].pos - sim.bodies[0].pos))
    for _ in range(200):
        sim.step()
    r1 = float(np.linalg.norm(sim.bodies[1].pos - sim.bodies[0].pos))
    assert r1 == pytest.approx(r0, rel=0.05)


def test_total_energy_regression() -> None:
    """Guards accidental changes to gravity / softening / integrator."""

    def factory() -> list[Body]:
        return [
            Body("Sun", 1.989e30, [0, 0, 0], [0, 0, 0], "yellow", 7e8, max_traj_points=200),
            Body(
                "Earth",
                5.97e24,
                [1.496e11, 0, 0],
                [0, 29780, 0],
                "blue",
                6.4e6,
                max_traj_points=200,
            ),
        ]

    sim = SolarSystemSim(
        86400.0,
        max_traj_points=200,
        softening_m=1e6,
        body_factory=factory,
    )
    for _ in range(20):
        sim.step()
    e = sim.compute_total_energy()
    assert e == pytest.approx(-2.647e33, rel=0.02)
