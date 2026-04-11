@echo off
cd /d "%~dp0testbed\src\Testbed.Visual\bin\Debug\net9.0"
Testbed.Visual.exe
if errorlevel 1 pause
