#ifndef SWE_HALOS_H
#define SWE_HALOS_H

#include <swe_config.h>

#ifdef USE_MPI

void exchange_SVM_halos(fType* H);

#endif

#endif
