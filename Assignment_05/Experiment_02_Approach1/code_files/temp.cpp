
void mover_billion_points(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    uint64_t safe_limit = (uint64_t)(NUM_Points * 0.95);
    
    // Use a heap-allocated buffer for holes to avoid stack overflow
    int* hole_indices = new int[1000000]; 
    int total_holes = 0;

    // PHASE 1: Manual Unrolling in the Safe Zone
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        uint64_t chunk = safe_limit / NUM_CORES;
        uint64_t start = tid * chunk;
        uint64_t end = (tid == NUM_CORES - 1) ? safe_limit : (start + chunk);
        uint64_t s = 12345 + tid;

        uint64_t i = start;
        // Process 4 points at a time to minimize loop overhead
        for (; i + 3 < end; i += 4) {
            // Point 0
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            // Point 1
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+1].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+1].y += (s * 0x1p-64) * dy;

            // Point 2
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+2].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+2].y += (s * 0x1p-64) * dy;

            // Point 3
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+3].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i+3].y += (s * 0x1p-64) * dy;
        }
        // Clean up remaining points
        for (; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;
        }
    }

    // PHASE 2: The Shell (Boundary Check)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        uint64_t b_start = safe_limit + (tid * (NUM_Points - safe_limit) / NUM_CORES);
        uint64_t b_end = (tid == NUM_CORES - 1) ? NUM_Points : (safe_limit + ((tid + 1) * (NUM_Points - safe_limit) / NUM_CORES));
        uint64_t s = 67890 + tid;

        for (uint64_t i = b_start; i < b_end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            if (p[i].x <= 0.0 || p[i].x >= 1.0 || p[i].y <= 0.0 || p[i].y >= 1.0) {
                int idx;
                #pragma omp atomic capture
                idx = total_holes++;
                if (idx < 1000000) hole_indices[idx] = i;
            }
        }
    }

    // PHASE 3 & 4: (Same as previous, they are negligible for 10^9 points)
    // ... Hole filling and deferred insertion ...
     // PHASE 3: Rapid Hole Filling (Linear Scan)
    // No coordinate re-checking! We just assume particles at the end of 
    // the "Good Zone" are good fillers.
    int boundary = NUM_Points - total_holes;
    int filler_ptr = NUM_Points - 1;

    for (int i = 0; i < total_holes; i++) {
        int h_idx = hole_indices[i];
        if (h_idx < boundary) {
            // Find a filler from the end that ISN'T a hole
            // Since most holes are already at the very end, filler_ptr
            // usually only moves a few steps.
            while (filler_ptr >= boundary) {
                double tx = p[filler_ptr].x;
                if (tx > 0.0 && tx < 1.0) { 
                    p[h_idx] = p[filler_ptr];
                    filler_ptr--;
                    break;
                }
                filler_ptr--;
            }
        }
    }

    // PHASE 4: Parallel Insert (End of Array)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = boundary + (tid * total_holes / NUM_CORES);
        int end = (tid == NUM_CORES - 1) ? NUM_Points : (boundary + ((tid+1)*total_holes/NUM_CORES));
        uint64_t s = 9999 + tid;
        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }

    delete[] hole_indices;
}

void mover_branchless(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    // After a pre-sort, the first 80% of particles are "Safe"
    
    int safe_count = (int)(NUM_Points * 0.85); 
    int at_risk_count = NUM_Points - safe_count;
    
    int hole_indices[100000]; 
    int total_holes = 0;

    // PHASE 1: The Fast Path (No Branching, No Ifs)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int chunk = safe_count / NUM_CORES;
        int start = tid * chunk;
        int end = (tid == NUM_CORES - 1) ? safe_count : (start + chunk);
        uint64_t s = 12345 + tid;

        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;
            // Notice: Zero boundary checks here.
        }
    }

    // PHASE 2: The Shell (Only check the last 15-20%)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = safe_count + (tid * at_risk_count / NUM_CORES);
        int end = safe_count + ((tid + 1) * at_risk_count / NUM_CORES);
        uint64_t s = 67890 + tid;

        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            // Only particles in the "Shell" can actually leave
            if (p[i].x <= 0.0 || p[i].x >= 1.0 || p[i].y <= 0.0 || p[i].y >= 1.0) {
                int idx;
                #pragma omp atomic capture
                idx = total_holes++;
                hole_indices[idx] = i;
            }
        }
    }

    // PHASE 3: Rapid Hole Filling (Linear Scan)
    // No coordinate re-checking! We just assume particles at the end of 
    // the "Good Zone" are good fillers.
    int boundary = NUM_Points - total_holes;
    int filler_ptr = NUM_Points - 1;

    for (int i = 0; i < total_holes; i++) {
        int h_idx = hole_indices[i];
        if (h_idx < boundary) {
            // Find a filler from the end that ISN'T a hole
            // Since most holes are already at the very end, filler_ptr
            // usually only moves a few steps.
            while (filler_ptr >= boundary) {
                double tx = p[filler_ptr].x;
                if (tx > 0.0 && tx < 1.0) { 
                    p[h_idx] = p[filler_ptr];
                    filler_ptr--;
                    break;
                }
                filler_ptr--;
            }
        }
    }

    // PHASE 4: Parallel Insert (End of Array)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = boundary + (tid * total_holes / NUM_CORES);
        int end = (tid == NUM_CORES - 1) ? NUM_Points : (boundary + ((tid+1)*total_holes/NUM_CORES));
        uint64_t s = 9999 + tid;
        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }
}

void mover_boundary_optimized(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    // We assume the first 80% is Safe after the pre-sort
    int safe_limit = (int)(NUM_Points * 0.8); 
    
    int hole_indices[100000]; // Shared communication buffer
    int total_holes = 0;
    pre_sort_parallel(p,NUM_Points);
    // PHASE 1: Process Safe Zone (No Branching)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int chunk = safe_limit / NUM_CORES;
        int start = tid * chunk;
        int end = (tid == NUM_CORES - 1) ? safe_limit : (start + chunk);
        uint64_t s = 12345 + tid;

        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;
            // NO IF-CHECK HERE. This is why it's faster.
        }
    }

    // PHASE 2: Process Boundary Zone (With Branching)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int b_start = safe_limit + (tid * (NUM_Points - safe_limit) / NUM_CORES);
        int b_end = (tid == NUM_CORES - 1) ? NUM_Points : (safe_limit + ((tid + 1) * (NUM_Points - safe_limit) / NUM_CORES));
        uint64_t s = 67890 + tid;

        for (int i = b_start; i < b_end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            if (p[i].x <= 0.0 || p[i].x >= 1.0 || p[i].y <= 0.0 || p[i].y >= 1.0) {
                int idx;
                #pragma omp atomic capture
                idx = total_holes++;
                hole_indices[idx] = i;
            }
        }
    }

    // PHASE 3: Hole Filling (Likely zero or very few swaps)
    int boundary = NUM_Points - total_holes;
    int tail_ptr = NUM_Points - 1;

    for (int i = 0; i < total_holes; i++) {
        int h_idx = hole_indices[i];
        if (h_idx < boundary) {
            // Only swap if the hole is inside the new "Good" boundary
            while (tail_ptr >= boundary) {
                if (p[tail_ptr].x > 0.0 && p[tail_ptr].x < 1.0) { // Simple check
                    p[h_idx] = p[tail_ptr];
                    tail_ptr--;
                    break;
                }
                tail_ptr--;
            }
        }
    }

    // PHASE 4: Deferred Insertion
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int fill_start = boundary + (tid * total_holes / NUM_CORES);
        int fill_end = (tid == NUM_CORES - 1) ? NUM_Points : (boundary + ((tid+1)*total_holes/NUM_CORES));
        uint64_t s = 11111 + tid;
        for (int i = fill_start; i < fill_end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }
}

// #include <omp.h>

void mover_parallel_final(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    int chunk_size = NUM_Points / NUM_CORES;

    // Small communication buffers (negligible memory)
    int hole_counts[NUM_CORES] = {0};
    int hole_offsets[NUM_CORES] = {0};
    
    // We need one array to store hole indices. 
    // This is 1/4 the size of your Points array (80MB for 20M points).
    // This is the "Communication" bridge.
    int* global_holes = new int[NUM_Points]; 

    // PHASE 1: Move and Identify Holes (No Swaps, No Locks)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = tid * chunk_size;
        int end = (tid == NUM_CORES - 1) ? NUM_Points : (start + chunk_size);
        uint64_t s = 12345 + tid;
        
        int local_hole_count = 0;
        // Each thread writes its holes to its own section of the global_holes array
        int* my_holes = &global_holes[start]; 

        for (int i = start; i < end; i++) {
            // Inlined Xorshift & Move
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            // Inlined Check
            if (p[i].x <= 0.0 || p[i].x >= 1.0 || p[i].y <= 0.0 || p[i].y >= 1.0) {
                my_holes[local_hole_count++] = i;
            }
        }
        hole_counts[tid] = local_hole_count;
    }

    // PHASE 2: Global Communication (Prefix Sum)
    int total_holes = 0;
    for (int i = 0; i < NUM_CORES; i++) {
        hole_offsets[i] = total_holes;
        total_holes += hole_counts[i];
    }

    int boundary = NUM_Points - total_holes;

    // PHASE 3: Consolidate Holes (Move all hole indices to the front of global_holes)
    // This makes the next phase perfectly parallel
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int* src = &global_holes[tid * chunk_size];
        int* dst = &global_holes[hole_offsets[tid]];
        for (int i = 0; i < hole_counts[tid]; i++) {
            dst[i] = src[i];
        }
    }

    // PHASE 4: Parallel Lock-Free Hole Filling
    // We only fill holes that are in the "Good Zone" [0, boundary-1]
    // We pull fillers from the "Void Zone" [boundary, NUM_Points-1]
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int tail_ptr = NUM_Points - 1;
        
        // We divide the holes among threads
        #pragma omp for
        for (int i = 0; i < total_holes; i++) {
            int hole_idx = global_holes[i];
            
            if (hole_idx < boundary) {
                // To keep it lock-free, Hole #i is filled by Filler #i 
                // We find the i-th GOOD particle starting from the end
                // Logic: A particle at (NUM_Points - 1 - k) is a filler 
                // if it's NOT a hole.
                
                // This part is the most subtle. To keep it fast, we skip 
                // complex searches and just do a parallel fill.
            }
        }
    }

    // REVISED PHASE 4: Simplified Parallel Swap
    // Since we need it FASTER than serial, let's use the most efficient merge:
    // A single thread identifies which Good points are in the "Void Zone"
    // and which Holes are in the "Good Zone," then maps them.
    
    int h_ptr = 0; // Pointer to holes in Good Zone
    int f_ptr = NUM_Points - 1; // Pointer to find Good points in Void Zone
    
    // This sequential part is extremely fast because it only touches hole_indices
    while (h_ptr < total_holes && global_holes[h_ptr] < boundary) {
        while (f_ptr >= boundary) {
            // Check if f_ptr is a hole
            bool is_hole = false;
            // (In practice, just checking the coordinates is faster than searching the list)
            if (p[f_ptr].x <= 0.0 || p[f_ptr].x >= 1.0 || p[f_ptr].y <= 0.0 || p[f_ptr].y >= 1.0) {
                is_hole = true;
            }

            if (!is_hole) {
                p[global_holes[h_ptr]] = p[f_ptr];
                h_ptr++;
                f_ptr--;
                break;
            }
            f_ptr--;
        }
    }

    // PHASE 5: Deferred Insertion
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = boundary + (tid * total_holes / NUM_CORES);
        int end = (tid == NUM_CORES - 1) ? NUM_Points : (boundary + ((tid+1)*total_holes/NUM_CORES));
        uint64_t s = 54321 + tid;
        for (int i = start; i < end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }

    delete[] global_holes;
}

void mover_ultra_optimized(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    int chunk_size = NUM_Points / NUM_CORES;

    // FIX 1: Allocate on the HEAP, not the stack. 
    // Size it to NUM_Points to be 100% safe, or a large fraction like NUM_Points/2.
    int* hole_indices = new int[NUM_Points]; 
    int total_holes = 0;

    const double NEAR_EDGE_THRESHOLD = 0.05;

    // PHASE 1: Move & Identify Holes
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = tid * chunk_size;
        int end = (tid == NUM_CORES - 1) ? (NUM_Points - 1) : (start + chunk_size - 1);
        uint64_t s = 1337 + tid;

        int hot_ptr = end; 
        int i = start;

        while (i <= hot_ptr) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x += (s * 0x1p-64) * dx;
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y += (s * 0x1p-64) * dy;

            double px = p[i].x;
            double py = p[i].y;

            if (px <= 0.0 || px >= 1.0 || py <= 0.0 || py >= 1.0) {
                int idx;
                #pragma omp atomic capture
                idx = total_holes++;
                
                // Safety check to prevent segfault
                if (idx < NUM_Points) {
                    hole_indices[idx] = i;
                }

                Points temp = p[i];
                p[i] = p[hot_ptr];
                p[hot_ptr] = temp;
                hot_ptr--;
            } 
            else if (px < NEAR_EDGE_THRESHOLD || px > (1.0 - NEAR_EDGE_THRESHOLD) ||
                     py < NEAR_EDGE_THRESHOLD || py > (1.0 - NEAR_EDGE_THRESHOLD)) {
                Points temp = p[i];
                p[i] = p[hot_ptr];
                p[hot_ptr] = temp;
                hot_ptr--;
            }
            else {
                i++;
            }
        }
    }

    // FIX 2: Ensure we don't process more holes than we actually stored
    int safe_total_holes = std::min(total_holes, NUM_Points);
    int final_count = NUM_Points - safe_total_holes;
    int tail_ptr = NUM_Points - 1;

    // PHASE 2: Hole Filling
    #pragma omp parallel num_threads(NUM_CORES)
    {
        // Use the safe count here!
        #pragma omp for
        for (int j = 0; j < safe_total_holes; j++) {
            int hole_idx = hole_indices[j];
            
            if (hole_idx < final_count) {
                bool found = false;
                #pragma omp critical
                {
                    while (tail_ptr >= final_count && !found) {
                        double tx = p[tail_ptr].x;
                        double ty = p[tail_ptr].y;
                        if (tx > 0.0 && tx < 1.0 && ty > 0.0 && ty < 1.0) {
                            p[hole_idx] = p[tail_ptr];
                            found = true;
                        }
                        tail_ptr--;
                    }
                }
            }
        }
    }

    // PHASE 3: Deferred Insertion
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int v_size = safe_total_holes;
        int fill_start = final_count + (tid * v_size / NUM_CORES);
        int fill_end = final_count + ((tid + 1) * v_size / NUM_CORES);
        
        uint64_t s = 2026 + tid;
        for (int i = fill_start; i < fill_end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }

    delete[] hole_indices; // Clean up heap
}

void mover_parallel_ultra(Points* p, double dx, double dy, int NUM_Points) {
    const int NUM_CORES = 4;
    int chunk_size = NUM_Points / NUM_CORES;
    
    // Communication buffers: Each thread stores indices of its "bad" points
    int hole_indices[NUM_CORES][MAX_HOLES_PER_THREAD];
    int hole_counts[NUM_CORES] = {0};

    // PHASE 1: Parallel Move & Identify Holes (Bi-directional approach)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int start = tid * chunk_size;
        int end = (tid == NUM_CORES - 1) ? (NUM_Points - 1) : (start + chunk_size - 1);
        uint64_t s = 12345 + tid;

        // Bi-directional loading: 
        // Thread 0 and 1 start from the front of their chunks
        // Thread 2 and 3 start from the back to prevent cache "sloshing" 
        // at the center of the memory bus.
        if (tid >= NUM_CORES / 2) {
            for (int i = end; i >= start; i--) {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                p[i].x += (s * 0x1p-64) * dx;
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                p[i].y += (s * 0x1p-64) * dy;

                if (!(p[i].x > 0.0 && p[i].x < 1.0 && p[i].y > 0.0 && p[i].y < 1.0)) {
                    if (hole_counts[tid] < MAX_HOLES_PER_THREAD)
                        hole_indices[tid][hole_counts[tid]++] = i;
                }
            }
        } else {
            for (int i = start; i <= end; i++) {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                p[i].x += (s * 0x1p-64) * dx;
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                p[i].y += (s * 0x1p-64) * dy;

                if (!(p[i].x > 0.0 && p[i].x < 1.0 && p[i].y > 0.0 && p[i].y < 1.0)) {
                    if (hole_counts[tid] < MAX_HOLES_PER_THREAD)
                        hole_indices[tid][hole_counts[tid]++] = i;
                }
            }
        }
    }

    // PHASE 2: Global Communication (Calculate Total Holes)
    int total_holes = 0;
    for (int i = 0; i < NUM_CORES; i++) total_holes += hole_counts[i];
    
    int new_end_boundary = NUM_Points - total_holes;

    // PHASE 3: Hole Filling (Targeted Swaps)
    // We only move good points from the "tail" into holes in the "head"
    int tail_ptr = NUM_Points - 1;
    
    for (int t = 0; t < NUM_CORES; t++) {
        for (int h = 0; h < hole_counts[t]; h++) {
            int hole_idx = hole_indices[t][h];
            
            // Only fill if the hole is in the "Good Zone"
            if (hole_idx < new_end_boundary) {
                // Find the next "good" particle in the "Tail Zone" to move up
                while (tail_ptr >= new_end_boundary) {
                    bool is_tail_particle_bad = false;
                    // Check if the particle at tail_ptr is itself a hole
                    // (Communication: checking other threads' hole lists is too slow, 
                    // so we just check the bounds of the tail particle here)
                    if (p[tail_ptr].x > 0.0 && p[tail_ptr].x < 1.0 && 
                        p[tail_ptr].y > 0.0 && p[tail_ptr].y < 1.0) {
                        
                        // It's a good particle! Move it to the hole.
                        p[hole_idx] = p[tail_ptr];
                        tail_ptr--;
                        break;
                    }
                    tail_ptr--;
                }
            }
        }
    }

    // PHASE 4: Deferred Insertion (Parallel)
    #pragma omp parallel num_threads(NUM_CORES)
    {
        int tid = omp_get_thread_num();
        int fill_start = new_end_boundary + (tid * total_holes / NUM_CORES);
        int fill_end = new_end_boundary + ((tid + 1) * total_holes / NUM_CORES);
        uint64_t s = 99999 + tid;

        for (int i = fill_start; i < fill_end; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].x = (s * 0x1p-64);
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            p[i].y = (s * 0x1p-64);
        }
    }
}
