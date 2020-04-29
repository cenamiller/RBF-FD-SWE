// Written by Samuel Elliott, Summer 2017. Last updated- Richelle Streater, Summer 2018

#include <profiling.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#else
#include <sys/time.h>
#endif

extern timing_struct local_timer;

extern int mpi_rank;
extern int mpi_size;

// Time measured in microseconds
double getTime(){
    
#ifdef USE_MPI
  return MPI_Wtime();
#else
  struct timeval time;
  int err = gettimeofday(&time, 0);
  return ((double) (time.tv_sec * 1000000 + time.tv_usec)) * 1.0e-6;
#endif
    
}

void init_timer(timing_struct *timer){
  timer->t_rhs = 0.0;
  timer->t_update = 0.0;
}
