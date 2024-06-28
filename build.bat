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

REM Run the executable
.\Release\MyProject.exe
