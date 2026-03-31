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

void interpolation(double *mesh_value, Points *points)
{
    double inv_dx = (double)NX;
    double inv_dy = (double)NY;

    for (int p = 0; p < NUM_Points; p++)
    {
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

        mesh_value[base_idx] += wx_m * wy_m;
        mesh_value[base_idx + 1] += lx * wy_m;
        mesh_value[base_idx + GRID_X] += wx_m * ly;
        mesh_value[base_idx + GRID_X + 1] += lx * ly;
    }
}

void mover_serial_1(Points *points, double deltaX, double deltaY)
{
    int valid_count = 0;
    uint32_t state = 123456789;

    for(int p = 0; p < NUM_Points; p++)
    {
        double x = points[p].x;
        double y = points[p].y;

        double rx = rng_uniform01(&state);
        double ry = rng_uniform01(&state);

        double newx = x + (2.0 * rx - 1.0) * deltaX;
        double newy = y + (2.0 * ry - 1.0) * deltaY;

        if(newx >= 0.0 && newx <= 1.0 && newy >= 0.0 && newy <= 1.0)
        {
            points[valid_count].x = newx;
            points[valid_count].y = newy;
            valid_count++;
        }
    }

    for(int p = valid_count; p < NUM_Points; p++)
    {
        points[p].x = rng_uniform01(&state);
        points[p].y = rng_uniform01(&state);
    }
}

void mover_serial_2(Points *points, double deltaX, double deltaY)
{
    uint32_t state = 987654321;

    for(int p = 0; p < NUM_Points; p++)
    {
        double x = points[p].x;
        double y = points[p].y;

        double rx = rng_uniform01(&state);
        double ry = rng_uniform01(&state);

        double newx = x + (2.0 * rx - 1.0) * deltaX;
        double newy = y + (2.0 * ry - 1.0) * deltaY;

        if(newx >= 0.0 && newx <= 1.0 && newy >= 0.0 && newy <= 1.0)
        {
            points[p].x = newx;
            points[p].y = newy;
        }
        else
        {
            points[p].x = rng_uniform01(&state);
            points[p].y = rng_uniform01(&state);
        }
    }
}

void mover_parallel_1(Points *points, double deltaX, double deltaY)
{
    const double RNG_NORM = 2.3283064365386963e-10;

    #pragma omp parallel for schedule(static)
    for(int p = 0; p < NUM_Points; p++)
    {
        uint32_t state = 314159265 ^ (omp_get_thread_num() * 1999) ^ p;

        double x = points[p].x;
        double y = points[p].y;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        double rx = (double)state * RNG_NORM;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        double ry = (double)state * RNG_NORM;

        double newx = x + (2.0 * rx - 1.0) * deltaX;
        double newy = y + (2.0 * ry - 1.0) * deltaY;

        if(newx >= 0.0 && newx <= 1.0 && newy >= 0.0 && newy <= 1.0)
        {
            points[p].x = newx;
            points[p].y = newy;
        }
        else
        {
            points[p].x = -1.0;
        }
    }

    int valid_count = 0;
    for(int p = 0; p < NUM_Points; p++)
    {
        if(points[p].x != -1.0)
        {
            if (valid_count != p) {
                points[valid_count].x = points[p].x;
                points[valid_count].y = points[p].y;
            }
            valid_count++;
        }
    }

    #pragma omp parallel for schedule(static)
    for(int p = valid_count; p < NUM_Points; p++)
    {
        uint32_t state = 271828182 ^ (omp_get_thread_num() * 2026) ^ p;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        points[p].x = (double)state * RNG_NORM;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        points[p].y = (double)state * RNG_NORM;
    }
}

void mover_parallel_2(Points *points, double deltaX, double deltaY)
{
    const double RNG_NORM = 2.3283064365386963e-10;

    #pragma omp parallel for schedule(static)
    for(int p = 0; p < NUM_Points; p++)
    {
        uint32_t state = 314159265 ^ (omp_get_thread_num() * 1999) ^ p;

        double x = points[p].x;
        double y = points[p].y;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        double rx = (double)state * RNG_NORM;

        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        double ry = (double)state * RNG_NORM;

        double newx = x + (2.0 * rx - 1.0) * deltaX;
        double newy = y + (2.0 * ry - 1.0) * deltaY;

        if(newx >= 0.0 && newx <= 1.0 && newy >= 0.0 && newy <= 1.0)
        {
            points[p].x = newx;
            points[p].y = newy;
        }
        else
        {
            state ^= state << 13; state ^= state >> 17; state ^= state << 5;
            points[p].x = (double)state * RNG_NORM;

            state ^= state << 13; state ^= state >> 17; state ^= state << 5;
            points[p].y = (double)state * RNG_NORM;
        }
    }
}

void save_mesh(double *mesh_value)
{
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd)
    {
        printf("Error creating Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < GRID_Y; i++)
    {
        for (int j = 0; j < GRID_X; j++)
        {
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}