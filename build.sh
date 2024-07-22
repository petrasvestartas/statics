#!/bin/bash
# Change directory to the directory containing this script
cd "$(dirname "$0")"

# Create a directory for the build files
mkdir -p build
cd build

# Configure the project. Adjust the generator if needed.
cmake ..

# Build the project with the specified configuration
cmake --build . --config Release

# Run the executable
./MyProject



#!/bin/bash

echo "Current working directory: $(pwd)"

# Use an absolute path for BASE_DIR if the script's execution location might vary
# Example: BASE_DIR="/absolute/path/to/your/project/build/examples/chapter2"
BASE_DIR="./examples/chapter2"

# Check if the base directory exists
if [ -d "$BASE_DIR" ]; then
    # Iterate over each subdirectory in the base directory
    for SUBDIR in "$BASE_DIR"/*; do
        # Check if the subdirectory is a directory
        if [ -d "$SUBDIR" ]; then
            # Iterate over each executable in the subdirectory
            for example in "$SUBDIR"/*; do
                # Check if the file is an executable
                if [ -x "$example" ]; then
                    # Print in blue color
                    echo -e "\033[34mRunning example: $example\033[0m"
                    # Run the example
                    "$example"
                else
                    echo "Skipping non-executable file: $example"
                fi
            done
        fi
    done
else
    echo "Base directory does not exist: $BASE_DIR"
fi