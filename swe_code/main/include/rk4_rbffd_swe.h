#ifndef RK4_RBFFD_SWE_H
#define RK4_RBFFD_SWE_H

#include <swe_config.h>

// Approximates one timestep of the constrained SWE PDE solution using RK4 timestepping and RBF-FD methods for the spatial differentiation

void RK_substep(fType* H, fType* K, fType* F, fType* D, int substep_id);

void copy_fp_arr(fType* dest, fType* src, int size);

void print_SVM(fType* H);

#endif
