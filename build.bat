@echo off
REM Change directory to the directory containing this script
cd %~dp0

REM Create a directory for the build files
if not exist build mkdir build
cd build

REM Configure the project using the Visual Studio 2022 generator
cmake -G "Visual Studio 17 2022" ..

REM Build the project with the specified configuration
cmake --build . --config Release

REM Run the main executable
.\Release\MyProject.exe

REM Assuming the executables are in the Release folder under each project directory in examples
set "BASE_DIR=%~dp0build\examples"

REM Check if the base directory exists
if exist "%BASE_DIR%" (
    echo Current working directory: %cd%
    REM Find all executable files in the base directory and its subdirectories
    for /r "%BASE_DIR%" %%X in (*.exe) do (
        echo Running example: %%X
        "%%X"
    )
) else (
    echo Base directory does not exist: %BASE_DIR%
)