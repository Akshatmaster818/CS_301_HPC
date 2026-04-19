#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h> 

#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    if (argc != 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }
                                                                                                                                                                                                  
    FILE *file = fopen(argv[1], "rb");
    if (!file) {
        printf("Error opening input file\n");
        exit(1);
    }

    fread(&NX, sizeof(int), 1, file);
    fread(&NY, sizeof(int), 1, file);

    fread(&NUM_Points, sizeof(int), 1, file);
    fread(&Maxiter, sizeof(int), 1, file);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;
    const long data_offset = 4L * sizeof(int);

    double *mesh_value = (double *) calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points = (Points *) calloc(NUM_Points, sizeof(Points));
    const int thread_counts[] = {2, 4, 8, 16};
    const int num_configs = sizeof(thread_counts) / sizeof(thread_counts[0]);

    omp_set_dynamic(0);

    for (int cfg = 0; cfg < num_configs; cfg++) {
        int num_threads = thread_counts[cfg];
        omp_set_num_threads(num_threads);

        if (fseek(file, data_offset, SEEK_SET) != 0) {
            printf("Error resetting input file\n");
            free(mesh_value);
            free(points);
            fclose(file);
            return 1;
        }

        memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));

        double **thread_meshes = (double **)malloc(num_threads * sizeof(double *));
        for (int t = 0; t < num_threads; t++) {
            thread_meshes[t] = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
        }

        printf("\n");
        printf("============================================================\n");
        printf("               PARALLEL INTERPOLATION BENCHMARK             \n");
        printf("============================================================\n");
        printf(" Configuration Details:\n");
        printf("  - Grid Dimensions : %d x %d\n", NX, NY);
        printf("  - Total Particles : %d\n", NUM_Points);
        printf("  - Max Iterations  : %d\n", Maxiter);
        printf("  - OpenMP Threads  : %d\n", num_threads);
        printf("------------------------------------------------------------\n");

        double total_time = 0.0;

        for (int iter = 0; iter < Maxiter; iter++) {

            read_points(file, points);

            double start = omp_get_wtime();

            interpolation_parallel(mesh_value, points, thread_meshes, num_threads);

            double end = omp_get_wtime();

            double iter_time = end - start;
            total_time += iter_time;

            printf("  [Iteration %2d/%2d]  Time: %.6lf seconds\n", iter + 1, Maxiter, iter_time);
        }

        printf("------------------------------------------------------------\n");
        printf("  Total Time        : %.6lf seconds\n", total_time);
        printf("  Average Iter Time : %.6lf seconds\n", total_time / Maxiter);
        printf("============================================================\n\n");

        if (cfg == num_configs - 1) {
            save_mesh(mesh_value);
        }

        for (int t = 0; t < num_threads; t++) {
            free(thread_meshes[t]);
        }
        free(thread_meshes);
    }

    free(mesh_value);
    free(points);
    fclose(file);

    return 0;
}