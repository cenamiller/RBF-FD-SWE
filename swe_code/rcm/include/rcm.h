#ifndef SWE_RCM
#define SWE_RCM

#include <swe_config.h>

// returns integer id array pointer to the rcm ordered mapping
int* rcm_mapping(GSMD_struct GSMD);

// prints info on the maximum bandwidth of the RCM reordered sparse DM matrix
void print_max_bandwidth(GSMD_struct GSMD);

void min_max(int* degree, int* min, int* max, const int Nnbr, const int offset);
void reorder_increasing(int* x, int* degree, const int Nnbr, const int Nnodes);
void find_degree(const int* idx, const int* mapped, int* degree, int* x, const int node_id, const int Nnbr);
void print_mapping(int* mapping, int Nnodes);

#endif
