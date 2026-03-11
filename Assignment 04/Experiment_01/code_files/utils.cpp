#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include "utils.h"

 static inline uint32_t xorshift32(uint32_t *state)
 {
     uint32_t x = *state;
     x ^= x << 13;
     x ^= x >> 17;
     x ^= x << 5;
     *state = x;
     return x;
 }

 static inline double rng_uniform01(uint32_t *state)
 {
     return (double)xorshift32(state) / (double)(4294967295.0);
 }

// Interpolation (Serial Code)
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

// Stochastic Mover (Serial Code) 
void mover_serial(Points *points, double deltaX, double deltaY) {
     for (int p = 0; p < NUM_Points; p++) {
         double x = points[p].x;
         double y = points[p].y;

         double newx, newy;
         do {
             double rx = (double)rand() / (double)RAND_MAX;
             double ry = (double)rand() / (double)RAND_MAX;
             double ddx = (2.0 * rx - 1.0) * deltaX;
             double ddy = (2.0 * ry - 1.0) * deltaY;
             newx = x + ddx;
             newy = y + ddy;
         } while (newx < 0.0 || newx > 1.0 || newy < 0.0 || newy > 1.0);

         points[p].x = newx;
         points[p].y = newy;
     }
}

// Stochastic Mover (Parallel Code) 
void mover_parallel(Points *points, double deltaX, double deltaY) {
 #pragma omp parallel
     {
         uint32_t seed = (uint32_t)time(NULL);
         seed ^= (uint32_t)(0x9E3779B9u * (uint32_t)(omp_get_thread_num() + 1));
         seed ^= (uint32_t)(uintptr_t)points;

 #pragma omp for schedule(static)
         for (int p = 0; p < NUM_Points; p++) {
             uint32_t local = seed ^ (uint32_t)(747796405u * (uint32_t)(p + 1));

             double x = points[p].x;
             double y = points[p].y;

             double newx, newy;
             do {
                 double rx = rng_uniform01(&local);
                 double ry = rng_uniform01(&local);
                 double ddx = (2.0 * rx - 1.0) * deltaX;
                 double ddy = (2.0 * ry - 1.0) * deltaY;
                 newx = x + ddx;
                 newy = y + ddy;
             } while (newx < 0.0 || newx > 1.0 || newy < 0.0 || newy > 1.0);

             points[p].x = newx;
             points[p].y = newy;
         }
     }
}

// Write mesh to file
void save_mesh(double *mesh_value) {

    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) {
        printf("Error creating Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}