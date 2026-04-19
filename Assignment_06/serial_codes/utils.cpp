#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include <omp.h>
// Serial interpolation 
void interpolation(double *mesh_value, Points *points) {
    double inv_dx = (double)NX; 
    double inv_dy = (double)NY;

    for (int p = 0; p < NUM_Points; p++) {
        double x = points[p].x;
        double y = points[p].y;

        int i = (int)(x * inv_dx);
        int j = (int)(y * inv_dy);
        i = (i >= NX) ? NX - 1 : i;
        j = (j >= NY) ? NY - 1 : j;

        double lx = x - (i * dx);
        double ly = y - (j * dy);

        double wx_m = dx - lx;
        double wy_m = dy - ly;

        int base_idx = j * GRID_X + i;

        mesh_value[base_idx]              += wx_m * wy_m; 
        mesh_value[base_idx + 1]          += lx * wy_m; 
        mesh_value[base_idx + GRID_X]     += wx_m * ly; 
        mesh_value[base_idx + GRID_X + 1] += lx * ly; 
    }
}

 
void interpolation_parallel(double *mesh_value, Points *points, double **thread_meshes, int num_threads) {
    const double l_dx = dx;
    const double l_dy = dy;
    const double l_inv_dx = (double)NX;
    const double l_inv_dy = (double)NY;
    const int l_NX = NX;
    const int l_NY = NY;
    const int l_GRID_X = GRID_X;
    
    const int num_grid_pts = GRID_X * GRID_Y;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double *local_mesh = thread_meshes[tid];

        // Reset the local mesh to 0 for this iteration very quickly
        memset(local_mesh, 0, num_grid_pts * sizeof(double));

        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; p++) {
            double x = points[p].x;
            double y = points[p].y;

            int i = (int)(x * l_inv_dx);
            int j = (int)(y * l_inv_dy);

            i = (i >= l_NX) ? l_NX - 1 : i;
            j = (j >= l_NY) ? l_NY - 1 : j;

            double lx = x - (i * l_dx);
            double ly = y - (j * l_dy);

            double wx_m = l_dx - lx;
            double wy_m = l_dy - ly;

            int base_idx = j * l_GRID_X + i;

            local_mesh[base_idx]                 += wx_m * wy_m;
            local_mesh[base_idx + 1]             += lx * wy_m;
            local_mesh[base_idx + l_GRID_X]      += wx_m * ly;
            local_mesh[base_idx + l_GRID_X + 1]  += lx * ly;
        }

        #pragma omp for schedule(static)
        for (int k = 0; k < num_grid_pts; k++) {
            double sum = 0.0;
            for (int t = 0; t < num_threads; t++) {
                sum += thread_meshes[t][k];
            }
            // Write the reduced sum into the global mesh
            mesh_value[k] += sum;
        }
    }
}


// Write mesh to file
void save_mesh(double *mesh_value) {
    FILE *fd = fopen("Mesh1.out", "w");
    if (!fd) {
        printf("Error creating Mesh1.out\n");
        exit(1);
    }

    // Note: ensure GRID_Y and GRID_X are accessible here. 
    // They are declared globally in main.cpp, so you need to declare them as extern 
    // at the top of utils.cpp if you haven't already.
    extern int GRID_X, GRID_Y;

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}
