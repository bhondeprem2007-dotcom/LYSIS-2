@echo off
cd /d "%~dp0"
echo ========================================
echo          LYSIS RESEARCH DESK
echo ========================================
py -m pip install -r requirements.txt
if errorlevel 1 (
  echo NumPy installation failed. Check your Python installation and internet connection.
  pause
  exit /b 1
)
py app.py
pause
