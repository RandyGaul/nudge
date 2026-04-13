@echo off
cd /d "%~dp0"

echo Configuring native DLLs...
cmake -S testbed\native -B testbed\native\build -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 goto :fail

echo Building native DLLs...
cmake --build testbed\native\build --config Release
if errorlevel 1 goto :fail

echo Building visual testbed...
dotnet build -c Release testbed\src\Testbed.Visual\Testbed.Visual.csproj -v q
if errorlevel 1 goto :fail

echo Launching...
cd testbed\src\Testbed.Visual\bin\Release\net9.0
Testbed.Visual.exe
if errorlevel 1 pause
goto :eof

:fail
echo Build failed.
pause
