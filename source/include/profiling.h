#ifndef SWE_PROFILING_H
#define SWE_PROFILING_H

#include <swe_config.h>

// Time measured in microseconds
double getTime();

timing_struct init_timer(int nattempts);

void process_profiling_results();

#endif
