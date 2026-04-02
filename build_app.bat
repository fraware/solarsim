@echo off
pip install pyinstaller
pyinstaller solarsim.spec
echo Build complete. Executable: dist\solarsim.exe
pause