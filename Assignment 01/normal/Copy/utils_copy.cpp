#include <math.h>
#include "utils_copy.h"

void vector_copy_operation(double *x, double *y, int Np) {
    // Copy: x[p] = y[p]
    for (int p = 0; p < Np; p++) {
        x[p] = y[p];

        // Prevent compiler from optimizing away the loop
        if (((double)p) == 333.333)
            dummy(p);
    }
}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
