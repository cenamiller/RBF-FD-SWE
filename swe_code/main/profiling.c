// Written by Samuel Elliott, Summer 2017. Last updated- Richelle Streater, Summer 2018
#ifdef CACHEQ
#include <cq.h>
#endif

#include <profiling.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#else
#include <time.h>
#endif

extern sim_params_struct sim_params;
extern timing_struct local_timer;

extern int mpi_rank;
extern int mpi_size;

// Time measured in microseconds
double getTime(){
    
#define NS_PER_SEC 1000000000L

#ifdef USE_MPI
    return MPI_Wtime();
#else
    struct timespec Time;
    clock_gettime(CLOCK_MONOTONIC, &Time);
    return ( ((double) (Time.tv_sec * NS_PER_SEC + Time.tv_nsec )) * 1E-9);
#endif
    
}

timing_struct init_timer(int nattempts) {

    timing_struct timer;
    timer.t_init = 0.0;
    timer.t_ocl = 0.0;
#ifdef CACHEQ
    timer.t_main = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
    timer.t_eval_rhs = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
    timer.t_mpi = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
    timer.t_eval_K = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
    timer.t_update_D = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
    timer.t_update_H = (double*) cq_calloc(TIMER_POOL, nattempts, sizeof(double));
#else
    timer.t_main = (double*) calloc(nattempts, sizeof(double));
    timer.t_eval_rhs = (double*) calloc(nattempts, sizeof(double));
    timer.t_mpi = (double*) calloc(nattempts, sizeof(double));
    timer.t_eval_K = (double*) calloc(nattempts, sizeof(double));
    timer.t_update_D = (double*) calloc(nattempts, sizeof(double));
    timer.t_update_H = (double*) calloc(nattempts, sizeof(double));
#endif
    return timer;
}


void process_profiling_results() {
	
    int nattempts = sim_params.nattempts;
    int nsteps = sim_params.nsteps;
    
    timing_struct global_sum_timer = init_timer(nattempts);
	
    #ifdef USE_MPI
    MPI_Reduce(&local_timer.t_init, &global_sum_timer.t_init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_main, global_sum_timer.t_main, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_eval_rhs, global_sum_timer.t_eval_rhs, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_eval_K, global_sum_timer.t_eval_K, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_update_D, global_sum_timer.t_update_D, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_update_H, global_sum_timer.t_update_H, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_timer.t_mpi, global_sum_timer.t_mpi, nattempts, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_timer.t_ocl, &global_sum_timer.t_ocl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
    #else
    global_sum_timer.t_init = local_timer.t_init;
    global_sum_timer.t_ocl = local_timer.t_ocl;
    memcpy(global_sum_timer.t_main, local_timer.t_main, sizeof(double) * nattempts);
    memcpy(global_sum_timer.t_eval_rhs, local_timer.t_eval_rhs, sizeof(double) * nattempts);
    memcpy(global_sum_timer.t_eval_K, local_timer.t_eval_K, sizeof(double) * nattempts);
    memcpy(global_sum_timer.t_update_D, local_timer.t_update_D, sizeof(double) * nattempts);
    memcpy(global_sum_timer.t_update_H, local_timer.t_update_H, sizeof(double) * nattempts);
    memcpy(global_sum_timer.t_mpi, local_timer.t_mpi, sizeof(double) * nattempts);
    #endif
    
    if (mpi_rank == 0) {
		
        timing_struct average_timer = init_timer(1);
        timing_struct min_timer = init_timer(1);
        timing_struct max_timer = init_timer(1);
		
        average_timer.t_init = global_sum_timer.t_init / mpi_size;
        average_timer.t_ocl = global_sum_timer.t_ocl / mpi_size;
        
        for (int i = 0; i < nattempts; i++) {

            double t_main_i = global_sum_timer.t_main[i] / (nsteps * mpi_size);
            average_timer.t_main[0] += t_main_i / nattempts;
            if (i == 0 || min_timer.t_main[0] > t_main_i) min_timer.t_main[0] = t_main_i;
            if (i == 0 || max_timer.t_main[0] < t_main_i) max_timer.t_main[0] = t_main_i;

            double t_eval_K_i = global_sum_timer.t_eval_K[i] / (nsteps * mpi_size);
            average_timer.t_eval_K[0] += t_eval_K_i / nattempts;
            if (i == 0 || min_timer.t_eval_K[0] > t_eval_K_i) min_timer.t_eval_K[0] = t_eval_K_i;
            if (i == 0 || max_timer.t_eval_K[0] < t_eval_K_i) max_timer.t_eval_K[0] = t_eval_K_i;

            double t_update_D_i = global_sum_timer.t_update_D[i] / (nsteps * mpi_size);
            average_timer.t_update_D[0] += t_update_D_i / nattempts;
            if (i == 0 || min_timer.t_update_D[0] > t_update_D_i) min_timer.t_update_D[0] = t_update_D_i;
            if (i == 0 || max_timer.t_update_D[0] < t_update_D_i) max_timer.t_update_D[0] = t_update_D_i;

            double t_update_H_i = global_sum_timer.t_update_H[i] / (nsteps * mpi_size);
            average_timer.t_update_H[0] += t_update_H_i / nattempts;
            if (i == 0 || min_timer.t_update_H[0] > t_update_H_i) min_timer.t_update_H[0] = t_update_H_i;
            if (i == 0 || max_timer.t_update_H[0] < t_update_H_i) max_timer.t_update_H[0] = t_update_H_i;
            
            double t_mpi_i = global_sum_timer.t_mpi[i] / (nsteps * mpi_size);
            average_timer.t_mpi[0] += t_mpi_i / nattempts;
            if (i == 0 || min_timer.t_mpi[0] > t_mpi_i) min_timer.t_mpi[0] = t_mpi_i;
            if (i == 0 || max_timer.t_mpi[0] < t_mpi_i) max_timer.t_mpi[0] = t_mpi_i;
        }

        if (nattempts > 5) {
        // Eliminate first two attempts from average for warmup time
           for (int i = 2; i < nattempts; i++) {
               double t_eval_rhs_i = global_sum_timer.t_eval_rhs[i] / (nsteps * mpi_size);
               average_timer.t_eval_rhs[0] += t_eval_rhs_i / (nattempts-2);
               if (i == 2 || min_timer.t_eval_rhs[0] > t_eval_rhs_i) min_timer.t_eval_rhs[0] = t_eval_rhs_i;
               if (i == 2 || max_timer.t_eval_rhs[0] < t_eval_rhs_i) max_timer.t_eval_rhs[0] = t_eval_rhs_i;
           }
        } else {
           for (int i = 0; i < nattempts; i++) {
               double t_eval_rhs_i = global_sum_timer.t_eval_rhs[i] / (nsteps * mpi_size);
               average_timer.t_eval_rhs[0] += t_eval_rhs_i / nattempts;
               if (i == 0 || min_timer.t_eval_rhs[0] > t_eval_rhs_i) min_timer.t_eval_rhs[0] = t_eval_rhs_i;
               if (i == 0 || max_timer.t_eval_rhs[0] < t_eval_rhs_i) max_timer.t_eval_rhs[0] = t_eval_rhs_i;
           }
	}

        timing_struct stddev2_timer = init_timer(1);

        for (int i = 0; i < nattempts; i++) {

            double t_main_i = global_sum_timer.t_main[i] / (nsteps * mpi_size);
            stddev2_timer.t_main[0] += pow(average_timer.t_main[0] - t_main_i, 2);

            double t_eval_rhs_i = global_sum_timer.t_eval_rhs[i] / (nsteps * mpi_size);
            stddev2_timer.t_eval_rhs[0] += pow(average_timer.t_eval_rhs[0] - t_eval_rhs_i, 2);

            double t_eval_K_i = global_sum_timer.t_eval_K[i] / (nsteps * mpi_size);
            stddev2_timer.t_eval_K[0] += pow(average_timer.t_eval_K[0] - t_eval_K_i, 2);

            double t_update_D_i = global_sum_timer.t_update_D[i] / (nsteps * mpi_size);
            stddev2_timer.t_update_D[0] += pow(average_timer.t_update_D[0] - t_update_D_i, 2);

            double t_update_H_i = global_sum_timer.t_update_H[i] / (nsteps * mpi_size);
            stddev2_timer.t_update_H[0] += pow(average_timer.t_update_H[0] - t_update_H_i, 2);

            double t_mpi_i = global_sum_timer.t_mpi[i] / (nsteps * mpi_size);
            stddev2_timer.t_mpi[0] += pow(average_timer.t_mpi[0] - t_mpi_i, 2);

        }

        printf("\n============================================================= Profiling Results ==============================================================\n\n");

        printf("Total Initialization Time (seconds): \t%e\n", average_timer.t_init);
        
        if (sim_params.OCL == 1) {       
            #ifdef USE_OCL
            printf("OpenCL Setup Time (seconds): \t%e\n", average_timer.t_ocl);
            #endif
        }
        printf("Main RK4 Loop (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_main[0], min_timer.t_main[0], max_timer.t_main[0], sqrt(stddev2_timer.t_main[0]) / nattempts);
        printf("Eval_Rhs      (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_eval_rhs[0], min_timer.t_eval_rhs[0], max_timer.t_eval_rhs[0], sqrt(stddev2_timer.t_eval_rhs[0]) / nattempts);
        printf("Eval_K        (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_eval_K[0], min_timer.t_eval_K[0], max_timer.t_eval_K[0], sqrt(stddev2_timer.t_eval_K[0]) / nattempts);
        printf("Update_D      (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_update_D[0], min_timer.t_update_D[0], max_timer.t_update_D[0], sqrt(stddev2_timer.t_update_D[0]) / nattempts);
        printf("Update_H      (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_update_H[0], min_timer.t_update_H[0], max_timer.t_update_H[0], sqrt(stddev2_timer.t_update_H[0]) / nattempts);
        
        #ifdef USE_MPI
        printf("MPI Overhead  (seconds/timestep) -> \tAverage: \t%e \tMin: \t%e \tMax: \t%e \tSTDDEV: \t%e\n", average_timer.t_mpi[0], min_timer.t_mpi[0], max_timer.t_mpi[0], sqrt(stddev2_timer.t_mpi[0]) / nattempts);
        #endif
        
        printf("\n==============================================================================================================================================\n");
    }
}
