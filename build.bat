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
    
    REM List directories first, then files, sorted alphabetically
    for /f "delims=" %%D in ('dir /b /ad "%BASE_DIR%" ^| sort') do (
        REM Find all executable files in each directory
        for /r "%BASE_DIR%\%%D" %%X in (*.exe) do (
            echo Running example: %%X
            "%%X"
        )
    )
    
    REM List executable files in the base directory, sorted alphabetically
    for /f "delims=" %%F in ('dir /b /a-d "%BASE_DIR%\*.exe" ^| sort') do (
        echo Running example: %BASE_DIR%\%%F
        "%BASE_DIR%\%%F"
    )
) else (
    echo Base directory does not exist: %BASE_DIR%
)