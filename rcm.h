#ifndef _RCM_
#define _RCM_

// returns integer id array pointer to the rcm ordered mapping
int* rcm_mapping(const int NNodes, const int NNbrs_Padded, const int NNbrs, int *idx);

// prints info on the maximum bandwidth of the RCM reordered sparse DM matrix
void rcm_print_max_bandwidth(const int NNodes, const int NNbrs_Padded, const int NNbrs, int *idx);

void rcm_min_max(int* degree, int* min, int* max, const int NNbrs, const int offset);
void rcm_reorder_increasing(int* x, int* degree, const int NNbrs, const int NNodes);
void rcm_find_degree(const int* idx, const int* mapped, int* degree, int* x, const int node_id, const int NNbrs_Padded, const int NNbrs);
void rcm_print_mapping(int* mapping, int NNodes);
void rcm_check_mapping(int* mapping, int NNodes);


#endif
