#ifndef PROFILING_H
#define PROFILING_H

typedef struct timing_struct {

  double t_rhs;
  double t_update;

} timing_struct;

// Time measured in microseconds

double getTime();
void init_timer(timing_struct* timer);
#endif
