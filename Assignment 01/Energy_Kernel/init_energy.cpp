#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "init_energy.h"

void init_vectors(int Np, double **x, double **y, double **v, double **E) {

    // Allocate memory for vectors
    *x = (double*) malloc(Np * sizeof(double));
    *y = (double*) malloc(Np * sizeof(double));
    *v = (double*) malloc(Np * sizeof(double));
    *E = (double*) malloc(Np * sizeof(double));

    for (int i = 0; i < Np; i++) {

        // Uniform random numbers in [0, 1] 
        (*x)[i] = (double) rand() / (double) RAND_MAX;
        (*y)[i] = (double) rand() / (double) RAND_MAX;
        (*v)[i] = (double) rand() / (double) RAND_MAX;

        (*E)[i] = 0.0;
    }
}

void free_vectors(double *x, double *y, double *v, double *E) {
    free(x);
    free(y);
    free(v);
    free(E);
}
