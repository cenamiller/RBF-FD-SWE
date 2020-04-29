#ifndef PROFILING_H
#define PROFILING_H

#include <stdint.h>
typedef struct timing_struct {

  double t_rhs;
  double t_update;

} timing_struct;

#ifdef CACHEQ
typedef struct CQ_timing_struct {
    
    uint64_t t_rhs;
    uint64_t t_update;
} CQ_timing_struct;
#endif


// Time measured in microseconds
#ifndef CACHEQ
double getTime();
#endif
void init_timer(timing_struct* timer);
#endif
