#pragma once

#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits> 

// Solves system of linear equations using Gaussian elimination
// Input matrix should be in augmented form [A|b] where A is the coefficient matrix
// and b is the right-hand side vector
// @param matrix - The augmented matrix [A|b] where the last column is b
// @return - Solution array of size N, where N is the number of variables
template<size_t N>
std::array<double, N> gaussian_elimination(std::array<std::array<double, N+1>, N>& matrix) {
    std::array<double, N> solution;

    // Forward elimination
    for (size_t i = 0; i < N; i++) {
        // Find pivot
        size_t pivot_row = i;
        double max_val = std::abs(matrix[i][i]);
        
        for (size_t j = i + 1; j < N; j++) {
            if (std::abs(matrix[j][i]) > max_val) {
                max_val = std::abs(matrix[j][i]);
                pivot_row = j;
            }
        }

        // Check if matrix is singular
        if (std::abs(max_val) < 1e-10) {
            // Return array filled with NaN to indicate no solution
            solution.fill(std::numeric_limits<double>::quiet_NaN());
            return solution;
        }

        // Swap rows if necessary
        if (pivot_row != i) {
            matrix[i].swap(matrix[pivot_row]);
        }

        // Eliminate column entries below pivot
        for (size_t j = i + 1; j < N; j++) {
            double factor = matrix[j][i] / matrix[i][i];
            for (size_t k = i; k < N + 1; k++) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        double sum = matrix[i][N];  // Get constant term
        for (size_t j = i + 1; j < N; j++) {
            sum -= matrix[i][j] * solution[j];
        }
        solution[i] = sum / matrix[i][i];
    }

    return solution;
}