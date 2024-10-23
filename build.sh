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

# Run the main executable
./MyProject

# Run all executables in the examples directory and its subdirectories
echo "Current working directory: $(pwd)"

# Use an absolute path for BASE_DIR if the script's execution location might vary
BASE_DIR="./examples"

# # Check if the base directory exists
# if [ -d "$BASE_DIR" ]; then
#     # Find all directories and files in the base directory and its subdirectories
#     find "$BASE_DIR" -type d -o -type f -executable | sort | while read -r item; do
#         if [ -d "$item" ]; then
#             # If the item is a directory, find all executables within it
#             find "$item" -type f -executable | sort | while read -r example; do
#                 # Print in blue color
#                 echo -e "\033[34mRunning example: $example\033[0m"
#                 # Run the example
#                 "$example"
#             done
#         elif [ -x "$item" ]; then
#             # If the item is an executable file, run it
#             echo -e "\033[34mRunning example: $item\033[0m"
#             "$item"
#         fi
#     done
# else
#     echo "Base directory does not exist: $BASE_DIR"
# fi