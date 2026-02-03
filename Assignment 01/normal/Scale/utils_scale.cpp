#include <math.h>
#include "utils_scale.h"

void vector_scale_operation(double *x, double *y, int Np) {
    // Scale: x[p] = a * y[p], where a = 2.5
    double a = 2.5;

    for (int p = 0; p < Np; p++) {
        x[p] = a * y[p];

        // Prevent compiler from optimizing away the loop
        if (((double)p) == 333.333)
            dummy(p);
    }
}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
