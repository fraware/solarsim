# Solar System Simulation

This project is a 3D simulation of the solar system that includes enhanced physics (adaptive time-stepping, advanced integrators, relativistic corrections, collision detection) and interactive controls using Tkinter and Matplotlib.

## Features

- Adaptive time-stepping and velocity Verlet integration.
- Realistic orbital parameters with a simplified relativistic correction for Mercury.
- Collision detection and merging of bodies.
- Interactive 3D visualization with camera controls.
- Timeline slider and data overlays (simulation time and total energy).

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/<your_username>/solar-system-sim.git
   cd solar-system-sim
   ```

Create a virtual environment and install dependencies:
bash
Copier
python -m venv venv
source venv/bin/activate # (venv\Scripts\activate on Windows)
pip install -r requirements.txt
Run the simulation:
bash
Copier
python enhanced_solar_system.py
Packaging as an Executable
(Optional) Use PyInstaller to create a standalone executable:

bash
Copier
pip install pyinstaller
pyinstaller --onefile enhanced_solar_system.py
The executable will be located in the dist/ folder.

License
This project is licensed under the MIT License. See the LICENSE file for details.
