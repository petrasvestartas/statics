# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

# Define the project and the languages used
project(MyProject VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable for the project. Include main.cpp and all source files from the src directory
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} main.cpp ${SOURCES})

# Include the src directory and its subdirectories to pick up the .hpp files
include_directories(${PROJECT_SOURCE_DIR}/src)

# If Point.hpp is in a subdirectory, explicitly include that directory as well
# include_directories(${PROJECT_SOURCE_DIR}/src/subdirectory)