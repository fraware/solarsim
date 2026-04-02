# PyInstaller spec for SolarSim (Tk + Matplotlib + Numba).
# Usage: pyinstaller solarsim.spec

import os

from PyInstaller.utils.hooks import collect_submodules

block_cipher = None
project_dir = os.path.dirname(os.path.abspath(SPEC))

src_dir = os.path.join(project_dir, "src")
entry = os.path.join(src_dir, "solarsim", "__main__.py")

hidden = ["tkinter", "matplotlib.backends.backend_tkagg"]
hidden += collect_submodules("matplotlib.backends")
hidden += collect_submodules("numba")

a = Analysis(
    [entry],
    pathex=[src_dir],
    binaries=[],
    datas=[],
    hiddenimports=hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name="solarsim",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
