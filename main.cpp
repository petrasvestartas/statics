#include <iostream>  // Includes the Standard Input Output Streams Library

#include "tests/test_arc.hpp"
#include "tests/test_intersection.hpp"
#include "tests/test_matrix.hpp"
#include "tests/test_pline.hpp"
#include "tests/test_point_and_vector.hpp"
#include "tests/test_vector.hpp"
#include "tests/test_xform.hpp"
#include "tests/test_gaussian_elimination.hpp"



int main() {

    std::cout << "Running tests..." << std::endl; // Outputs "Hello, World!" to the console
    test_vector_main();
    test_point_and_vector_main();
    test_arc_main();
    test_matrix_main();
    test_xform_main();
    test_pline_main(); // dont work
    test_intersection_main(); // not implemented
    test_gaussian_elimination_main();


    return 0;
}
