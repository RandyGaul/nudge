@echo off
REM Run tests with code coverage on bvh.c via OpenCppCoverage.
REM Generates HTML report in coverage/ and cobertura XML.

"C:\Program Files\OpenCppCoverage\OpenCppCoverage.exe" ^
  --sources src\bvh.c ^
  --export_type html:coverage ^
  --export_type cobertura:coverage.xml ^
  -- "%~dp0build\Debug\nudge_tests.exe"

echo.
for /f "tokens=2 delims==" %%a in ('findstr "line-rate" coverage.xml ^| findstr "coverage line"') do (
    echo BVH line coverage: %%a
)
echo Report: coverage\index.html
