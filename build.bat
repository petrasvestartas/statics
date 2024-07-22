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

REM Assuming the executables are in the Release folder under each project directory in examples\chapter2
set "BASE_DIR=%~dp0build\examples\chapter2"

REM Check if the base directory exists
if exist "%BASE_DIR%" (
    cd "%BASE_DIR%"
    for /d %%D in (*) do (
        cd "%%D"
        if exist "Release" (
            cd Release
            for %%X in (*.exe) do (
                if exist "%%X" (
                    echo Running example: %%X
                    .\%%X
                )
            )
            cd ..
        ) else (
            echo Release directory does not exist in: %%D
        )
        cd ..
    )
) else (
    echo Base directory does not exist: %BASE_DIR%
)