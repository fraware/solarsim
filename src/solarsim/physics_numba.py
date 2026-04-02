"""JIT-compiled N-body acceleration kernel (Plummer softening, optional Mercury tweak)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
from numba import njit

if TYPE_CHECKING:
    from solarsim.body import Body


@njit(cache=True)
def accelerations_kernel(
    pos: npt.NDArray[np.float64],
    vel: npt.NDArray[np.float64],
    mass: npt.NDArray[np.float64],
    mercury_idx: int,
    g: float,
    c_light: float,
    eps2: float,
) -> npt.NDArray[np.float64]:
    """Pairwise accelerations with softening ``sqrt(r^2 + eps^2)``.

    For body index ``mercury_idx`` (>= 0), applies an *illustrative* extra factor
    motivated by the GR magnitude of perihelion precession (not a full PN
    expansion). Documented on ``SolarSystemSim``.
    """
    n = pos.shape[0]
    acc = np.zeros((n, 3), dtype=np.float64)
    c2 = c_light * c_light
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            r2 = dx * dx + dy * dy + dz * dz + eps2
            inv_r = 1.0 / np.sqrt(r2)
            inv_r3 = inv_r / r2
            scale = g * mass[j] * inv_r3
            ax = scale * dx
            ay = scale * dy
            az = scale * dz
            if mercury_idx >= 0 and i == mercury_idx:
                lx = pos[i, 1] * vel[i, 2] - pos[i, 2] * vel[i, 1]
                ly = pos[i, 2] * vel[i, 0] - pos[i, 0] * vel[i, 2]
                lz = pos[i, 0] * vel[i, 1] - pos[i, 1] * vel[i, 0]
                l2 = lx * lx + ly * ly + lz * lz
                r2_geom = dx * dx + dy * dy + dz * dz
                if r2_geom > 0.0:
                    rel_corr = 1.0 + 3.0 * l2 / (r2_geom * c2)
                    ax *= rel_corr
                    ay *= rel_corr
                    az *= rel_corr
            acc[i, 0] += ax
            acc[i, 1] += ay
            acc[i, 2] += az
    return acc


def stack_state(bodies: list[Body]) -> tuple[npt.NDArray[np.float64], ...]:
    """Pack ``Body`` list into contiguous float64 arrays."""
    n = len(bodies)
    pos = np.zeros((n, 3), dtype=np.float64)
    vel = np.zeros((n, 3), dtype=np.float64)
    mass = np.zeros(n, dtype=np.float64)
    for i, b in enumerate(bodies):
        pos[i] = b.pos
        vel[i] = b.vel
        mass[i] = b.mass
    return pos, vel, mass
