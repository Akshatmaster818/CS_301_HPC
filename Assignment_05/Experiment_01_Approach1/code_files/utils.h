#define UTILS_H
#include <time.h>
#include "init.h"

void interpolation(double *mesh_value, Points *points);
void mover_serial(Points *points, double deltaX, double deltaY);
void mover_serial_del(Points* p,double dx,double dy);
void mover_serial_1(Points *points, double deltaX, double deltaY);
void mover_serial_2(Points *points, double deltaX, double deltaY);
void mover_parallel_1(Points *points, double deltaX, double deltaY);
void mover_parallel_2(Points *points, double deltaX, double deltaY);
// void mover_serial_del(Points* p,double dx,double dy);
// void mover_parallel_del(Points* p,double dx,double dy);
// // void mover_parallel_optimized(Points* P,double dx,double dy);
// void mover_no_opt_parallel(Points* p, double dx, double dy, int NUM_Points);
// void mover_parallel_optimized(Points* p, double dx, double dy, int NUM_Points);
// void mover_parallel_inplace(Points* p, double dx, double dy, int NUM_Points);
// void mover_parallel_ultra(Points* p, double dx, double dy, int NUM_Points);
// void mover_spatial_binning(Points* p, double dx, double dy, int NUM_Points);
// void mover_ultra_optimized(Points* p, double dx, double dy, int NUM_Points);
// void mover_ultimate_immediate(Points* p, double dx, double dy, int NUM_Points);
// void mover_parallel_final(Points* p, double dx, double dy, int NUM_Points);
// void mover_boundary_optimized(Points* p, double dx, double dy, int NUM_Points);
// void mover_branchless(Points* p, double dx, double dy, int NUM_Points);
// void pre_sort_parallel(Points* p, int NUM_Points);
// void mover_billion_points(Points* p, double dx, double dy, int NUM_Points);
// void mover_billion_nop(Points* p, double dx, double dy, int NUM_Points);
// void mover_parallel(Points *points, double deltaX, double deltaY);
// void mover_monolithic(Points* __restrict p, double dx, double dy, int NUM_Points);
// void mover_prefix_sum(Points* __restrict p, double dx, double dy, int NUM_Points);
// void mover_fenwick(Points* __restrict p, double dx, double dy, int NUM_Points);
void save_mesh(double *mesh_value);

#endif
