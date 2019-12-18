#ifndef RK4_RBFFD_SWE_H
#define RK4_RBFFD_SWE_H

#include <swe_config.h>
#ifndef CQ_POOL
#define CQ_POOL(x)
#endif

#define CQ_POOL_DF 3
#define CQ_POOL_HK 8

#define CQ_POOL_K2 7
#ifdef MOREPOOLS
#define CQ_POOL_K3 9
#define CQ_POOL_K4 10
#else
#define CQ_POOL_K3 CQ_POOL_HK
#define CQ_POOL_K4 CQ_POOL_K2
#endif

// Approximates one timestep of the constrained SWE PDE solution using RK4 timestepping and RBF-FD methods for the spatial differentiation
#ifdef CACHEQ

void RK_substep(fType CQ_POOL(CQ_POOL_HK) * H, fType CQ_POOL(CQ_POOL_HK) * K, 
                fType CQ_POOL(CQ_POOL_DF) * F, fType CQ_POOL(CQ_POOL_DF) * D, 
                int substep_id, PSMD_struct* LPSMD,
                timing_struct* local_timer);

void copy_fp_arr(fType CQ_POOL(CQ_POOL_DF) * dest, fType CQ_POOL(CQ_POOL_DF) * src, int size);

void print_SVM(fType CQ_POOL(CQ_POOL_HK) * H, PSMD_struct* LPSMD);

#else

void RK_substep(fType* H, fType* K, fType* F, fType* D, int substep_id);

void copy_fp_arr(fType* dest, fType* src, int size);

void print_SVM(fType* H);

#endif
#endif
