"""SolarSim: 3D N-body solar system simulation."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("solarsim")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = ["__version__"]
