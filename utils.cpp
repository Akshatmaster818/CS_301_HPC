#include <math.h>
#include "utils.h"

void vector_triad_operation(double *x, double *y, double *v, double *S, int Np) {

    for (int p = 0; p < Np; p++) {
        S[p] = x[p] + v[p] * y[p];

        // Prevent compiler from optimizing away the loop
        if (((double)p) == 333.333)
            dummy(p);

    }
}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
