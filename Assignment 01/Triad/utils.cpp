#include <math.h>
#include "utils.h"

// void vector_copy_operation(double *x, double *y, int Np) {}

// void vector_scale_operation(double *x, double *y, int Np) {}

// void vector_sum_operation(double *x, double *y, double *S, int Np) {}

void vector_triad_operation(double *x, double *y, double *v, double *S, int Np) {

    for (int p = 0; p < Np; p++) {
        S[p] = x[p] + v[p] * y[p];

         
        if (((double)p) == 333.333)
            dummy(p);

    }
}

// void energy_kernel_operation(double *v, double *E, int Np) {}

void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
