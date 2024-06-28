#include <iostream> // Includes the Standard Input Output Streams Library
#include "tests/test_vector.hpp"
#include "tests/test_point_and_vector.hpp"
#include "tests/test_arc.hpp"

// Main function - execution starts here
int main() {
    // std::cout << "Hello, World!" << std::endl; // Outputs "Hello, World!" to the console
    test_vector_main();
    test_point_and_vector_main();
    test_arc_main();
    return 0; // Returns 0 to signal the end of the program
}
