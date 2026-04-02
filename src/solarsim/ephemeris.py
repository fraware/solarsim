"""Initial conditions: built-in approximate elements and optional AstroPy ephemeris."""

from __future__ import annotations

import numpy as np

from solarsim.body import Body


def create_default_bodies(*, max_traj_points: int | None = 100_000) -> list[Body]:
    """Approximate masses, positions, and speeds (educational, not JPL precision)."""
    bodies: list[Body] = []
    kw = {"max_traj_points": max_traj_points}

    bodies.append(Body("Sun", 1.989e30, [0, 0, 0], [0, 0, 0], "yellow", 7e8, **kw))
    bodies.append(
        Body(
            "Mercury",
            3.3e23,
            [5.79e10, 0, 0],
            [0, 47400 * np.cos(np.deg2rad(7)), 47400 * np.sin(np.deg2rad(7))],
            "gray",
            2.4e6,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Venus",
            4.87e24,
            [1.082e11, 0, 0],
            [0, 35000 * np.cos(np.deg2rad(3.4)), 35000 * np.sin(np.deg2rad(3.4))],
            "orange",
            6e6,
            **kw,
        )
    )
    bodies.append(Body("Earth", 5.97e24, [1.496e11, 0, 0], [0, 29780, 0], "blue", 6.4e6, **kw))
    moon_dist = 3.84e8
    bodies.append(
        Body(
            "Moon",
            7.35e22,
            [1.496e11 + moon_dist, 0, 0],
            [0, 29780 + 1022, 0],
            "lightgray",
            1.7e6,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Mars",
            6.39e23,
            [2.279e11, 0, 0],
            [0, 24077 * np.cos(np.deg2rad(1.85)), 24077 * np.sin(np.deg2rad(1.85))],
            "red",
            3.4e6,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Jupiter",
            1.898e27,
            [7.78e11, 0, 0],
            [0, 13070 * np.cos(np.deg2rad(1.3)), 13070 * np.sin(np.deg2rad(1.3))],
            "saddlebrown",
            7e7,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Saturn",
            5.683e26,
            [1.43e12, 0, 0],
            [0, 9690 * np.cos(np.deg2rad(2.5)), 9690 * np.sin(np.deg2rad(2.5))],
            "goldenrod",
            6e7,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Uranus",
            8.681e25,
            [2.87e12, 0, 0],
            [0, 6810 * np.cos(np.deg2rad(0.77)), 6810 * np.sin(np.deg2rad(0.77))],
            "lightblue",
            2.5e7,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Neptune",
            1.024e26,
            [4.5e12, 0, 0],
            [0, 5430 * np.cos(np.deg2rad(1.77)), 5430 * np.sin(np.deg2rad(1.77))],
            "darkblue",
            2.5e7,
            **kw,
        )
    )
    bodies.append(
        Body(
            "Pluto",
            1.309e22,
            [5.9e12, 0, 0],
            [0, 4740 * np.cos(np.deg2rad(17)), 4740 * np.sin(np.deg2rad(17))],
            "purple",
            1.2e6,
            **kw,
        )
    )
    return bodies


def create_bodies_astropy(
    epoch: str = "2020-01-01",
    *,
    max_traj_points: int | None = 100_000,
) -> list[Body]:
    """Initial positions/velocities from AstroPy built-in solar-system ephemeris (SI units).

    Requires ``pip install 'solarsim[astropy]'``. Masses and display radii remain
    approximate constants; only kinematic state comes from ephemeris.
    """
    try:
        import astropy.units as u  # type: ignore[import-not-found]
        from astropy.coordinates import (  # type: ignore[import-not-found]
            get_body_barycentric_posvel,
            solar_system_ephemeris,
        )
        from astropy.time import Time  # type: ignore[import-not-found]
    except ImportError as exc:
        raise ImportError("Install optional dependency: pip install 'solarsim[astropy]'") from exc

    # kg; collision radii [m]; matplotlib colors
    catalog: dict[str, tuple[float, float, str]] = {
        "sun": (1.989e30, 7e8, "yellow"),
        "mercury": (3.3e23, 2.4e6, "gray"),
        "venus": (4.87e24, 6e6, "orange"),
        "earth": (5.97e24, 6.4e6, "blue"),
        "moon": (7.35e22, 1.7e6, "lightgray"),
        "mars": (6.39e23, 3.4e6, "red"),
        "jupiter": (1.898e27, 7e7, "saddlebrown"),
        "saturn": (5.683e26, 6e7, "goldenrod"),
        "uranus": (8.681e25, 2.5e7, "lightblue"),
        "neptune": (1.024e26, 2.5e7, "darkblue"),
        "pluto": (1.309e22, 1.2e6, "purple"),
    }

    display_names = {
        "sun": "Sun",
        "mercury": "Mercury",
        "venus": "Venus",
        "earth": "Earth",
        "moon": "Moon",
        "mars": "Mars",
        "jupiter": "Jupiter",
        "saturn": "Saturn",
        "uranus": "Uranus",
        "neptune": "Neptune",
        "pluto": "Pluto",
    }

    t = Time(epoch)
    bodies: list[Body] = []
    kw = {"max_traj_points": max_traj_points}

    with solar_system_ephemeris.set("builtin"):
        for key, (mass_kg, radius_m, color) in catalog.items():
            pos_vel = get_body_barycentric_posvel(key, t)
            r = pos_vel[0].get_xyz().to_value(u.m)
            v = pos_vel[1].get_xyz().to_value(u.m / u.s)
            name = display_names[key]
            bodies.append(Body(name, mass_kg, r, v, color, radius_m, **kw))

    return bodies
