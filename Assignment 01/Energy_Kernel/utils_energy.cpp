#include <math.h>
#include "utils_energy.h"

void energy_kernel_operation(double *v, double *E, int Np) {
    // Energy Kernel: E[p] = 0.5 * m * v[p]^2, where m = 1.0
    double m = 1.0;

    for (int p = 0; p < Np; p++) {
        E[p] = 0.5 * m * v[p] * v[p];

        // Prevent compiler from optimizing away the loop
        if (((double)p) == 333.333)
            dummy(p);
    }
}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
