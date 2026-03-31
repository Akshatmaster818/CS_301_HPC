#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "init.h"
#include "utils.h"

int pointsxy[3][2] = {{250, 100}, {500, 200}, {1000, 400}};
// int total_points[5] = {100, 10000, 1000000, 100000000, 1000000000};
// int total_points[1] = {100000000};
int total_points[1] = {14000000};

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv)
{

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            NX = pointsxy[i][0];
            NY = pointsxy[i][1];
            Maxiter = 10;
            NUM_Points = total_points[j];
            printf("Grid: %d x %d, Total Points: %d\n", NX, NY, NUM_Points);
            GRID_X = NX + 1;
            GRID_Y = NY + 1;
            dx = 1.0 / NX;
            dy = 1.0 / NY;

            omp_set_dynamic(0);

            int thread_counts[] = {1};
            int num_tests = sizeof(thread_counts) / sizeof(thread_counts[0]);

            double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
            Points *points = (Points *)calloc(NUM_Points, sizeof(Points));

            for (int t = 0; t < num_tests; t++)
            {
                int num_threads = thread_counts[t];
                omp_set_num_threads(num_threads);

                memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));
                initializepoints(points);

                printf("\nThreads: %d\n", num_threads);
                printf("Iter\tInterp\t\tMover\t\tTotal\n");

                for (int iter = 0; iter < Maxiter; iter++)
                {

                    double start_interp = omp_get_wtime();
                    interpolation(mesh_value, points);
                    double end_interp = omp_get_wtime();

                    double start_move = omp_get_wtime();
                    mover_serial_1(points, dx, dy);
                    double end_move = omp_get_wtime();

                    double interp_time = end_interp - start_interp;
                    double move_time = end_move - start_move;
                    double total = interp_time + move_time;

                    printf("%d\t%lf\t%lf\t%lf\n", iter + 1, interp_time, move_time, total);
                }
            }

            save_mesh(mesh_value);
            free(mesh_value);
            free(points);
        }
    }

    return 0;
}