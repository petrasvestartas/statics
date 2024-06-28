# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

# Define the project and the languages used
project(MyProject VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable for the project
add_executable(${PROJECT_NAME} main.cpp)

# Include the src directory to pick up the .hpp files
include_directories(src)

# Add all the source files from the src directory
file(GLOB SOURCES "src/*.cpp")
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
