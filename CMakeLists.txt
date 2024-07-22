# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

# Define the project and the languages used
project(MyProject VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable for the project. Include main.cpp and all source files from the src directory
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Create a static library from the sources instead of adding them directly to each executable
add_library(myproject_lib STATIC ${SOURCES})

# Add executable for the main project
add_executable(${PROJECT_NAME} main.cpp)
# Link the static library to the main executable
target_link_libraries(${PROJECT_NAME} myproject_lib)

# Include the src directory and its subdirectories to pick up the .hpp files
target_include_directories(myproject_lib PUBLIC ${PROJECT_SOURCE_DIR}/src)

# Create executables for each file in the examples/chapter2 directory
file(GLOB CHAPTER2_EXAMPLES "examples/chapter2/*.cpp")
foreach(EXAMPLE_FILE ${CHAPTER2_EXAMPLES})
    # Get the filename without the directory and extension
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
    # Create an executable for each example
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE})
    # Link the static library to the example executable
    target_link_libraries(${EXAMPLE_NAME} myproject_lib)
    # Set the output directory for the executables
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/chapter2/${EXAMPLE_NAME})
endforeach()