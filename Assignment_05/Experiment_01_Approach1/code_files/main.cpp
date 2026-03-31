#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

// Global variables
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int pointsxy[3][2] = {{250, 100}, {500, 200}, {1000, 400}};
int total_points[4] = {100, 10000, 1000000, 100000000};
// int total_points[1] = {100000000};
// int total_points[1] = {14000000};

int main(int argc, char **argv)
{

    for (int i = 0; i < 3; i++)
    {

        for (int j = 0; j < 4; j++)
        {
            // Fixed Parameters

            NX = pointsxy[i][0];
            NY = pointsxy[i][1];
            Maxiter = 10;
            NUM_Points = total_points[j];
            printf("Grid: %d x %d, Total Points: %d\n", NX, NY, NUM_Points);
            // Since Number of points will be 1 more than number of cells
            GRID_X = NX + 1;
            GRID_Y = NY + 1;
            dx = 1.0 / NX;
            dy = 1.0 / NY;

            // Fix Number of Threads
            // omp_set_num_threads(4);

            // Allocate memory for grid and Points
            double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
            Points *points = (Points *)calloc(NUM_Points, sizeof(Points));

            initializepoints(points);
            // pre_sort_parallel(points,NUM_Points);
            printf("Iter\tInterp\t\tMover\t\tTotal\n");
            for (int iter = 0; iter < Maxiter; iter++)
            {

                // Interpolation timing
                clock_t start_interp = clock();
                interpolation(mesh_value, points);
                clock_t end_interp = clock();
                // Mover timing
                clock_t start_move = clock();
                mover_serial_1(points, dx, dy);
                // mover_serial_del(points,dx,dy);
                // mover_parallel_del(points,dx,dy);
                // mover_ultimate_immediate(points,dx,dy,NUM_Points);
                // mover_serial_2(points,dx,dy);
                // mover_serial_del(points,dx,dy);
                clock_t end_move = clock();

                double interp_time = (double)(end_interp - start_interp) / CLOCKS_PER_SEC;
                double move_time = (double)(end_move - start_move) / CLOCKS_PER_SEC;
                double total = interp_time + move_time;

                printf("%d\t%lf\t%lf\t%lf\n", iter + 1, interp_time, move_time, total);
            }

            // Free memory
            save_mesh(mesh_value);
            free(mesh_value);
            free(points);
        }
    }

    return 0;
}