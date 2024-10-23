#include <iostream> // Includes the Standard Input Output Streams Library
#include "tests/test_vector.hpp"
#include "tests/test_point_and_vector.hpp"
#include "tests/test_arc.hpp"
#include "tests/test_limit_analysis.hpp"
#include "tests/test_matrix.hpp"
#include "tests/test_xform.hpp"
#include "tests/test_pline.hpp"

// Main function - execution starts here
int main() {
    // std::cout << "Hello, World!" << std::endl; // Outputs "Hello, World!" to the console
    test_vector_main();
    test_point_and_vector_main();
    test_arc_main();
    test_limit_analysis_main();
    test_matrix_main();
    test_xform_main();
    test_pline_main();
    return 0; // Returns 0 to signal the end of the program
}
