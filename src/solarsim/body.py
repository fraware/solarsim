from __future__ import annotations

from collections import deque

import numpy as np
import numpy.typing as npt


class Body:
    """Celestial body with state and bounded trajectory history."""

    def __init__(
        self,
        name: str,
        mass: float,
        pos: npt.ArrayLike,
        vel: npt.ArrayLike,
        color: str,
        radius: float,
        *,
        max_traj_points: int | None = 100_000,
    ) -> None:
        self.name = name
        self.mass = float(mass)
        self.pos = np.asarray(pos, dtype=np.float64)
        self.vel = np.asarray(vel, dtype=np.float64)
        self.color = color
        self.radius = float(radius)
        if max_traj_points is None:
            self.traj: list[npt.NDArray[np.float64]] | deque[npt.NDArray[np.float64]] = [
                self.pos.copy()
            ]
        else:
            self.traj = deque([self.pos.copy()], maxlen=max_traj_points)

    def append_traj(self, p: npt.NDArray[np.float64]) -> None:
        if isinstance(self.traj, deque):
            self.traj.append(p.copy())
        else:
            self.traj.append(p.copy())

    def extend_traj_from(self, other: Body) -> None:
        if isinstance(self.traj, deque) and isinstance(other.traj, deque):
            self.traj.extend(other.traj)
        elif isinstance(self.traj, deque):
            self.traj.extend(list(other.traj))
        else:
            self.traj.extend(list(other.traj))
