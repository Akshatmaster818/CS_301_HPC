#include <math.h>
#include "utils_add.h"

void vector_sum_operation(double *x, double *y, double *S, int Np) {
    // Add: S[p] = x[p] + y[p]
    for (int p = 0; p < Np; p++) {
        S[p] = x[p] + y[p];

        // Prevent compiler from optimizing away the loop
        if (((double)p) == 333.333)
            dummy(p);
    }
}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
