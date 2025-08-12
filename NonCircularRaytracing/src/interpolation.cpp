//Can you write a three point interpolation function in C++? --- IGNORE ---
#include "def.h"
long double three_point_interpolation(long double x, long double x0, long double x1, long double x2, long double f0, long double f1, long double f2) {
    // Calculate the coefficients for the Lagrange polynomial
    long double a = (f0 * (x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2));
    long double b = (f1 * (x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2));
    long double c = (f2 * (x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));

    // Return the interpolated value
    return a + b + c;
}
// This function uses the three-point interpolation formula to estimate the value at x based on the values
// at x0, x1, and x2 with their corresponding function values f0, f1, and f2.
// It is useful for estimating values in a smooth function where you have three known points.   