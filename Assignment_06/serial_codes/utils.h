#ifndef UTILS_H
#define UTILS_H
#include <time.h>
#include "init.h"

// PIC operations
void interpolation(double *mesh_value, Points *points);
void interpolation_parallel(double *mesh_value, Points *points, double **thread_meshes, int num_threads);
void parallel_interpolation_13(double* mesh_value, Points* points);
void parallel_interpolation_12(double* mesh_value, Points* points);
void save_mesh(double *mesh_value);

#endif
