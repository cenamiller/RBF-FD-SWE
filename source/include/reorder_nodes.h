#ifndef SWE_REORDER_NODES_H
#define SWE_REORDER_NODES_H

#include <swe_config.h>


// reorders any two dimensional floating point array to the new mapping (1D mapping must correspond to slowest growing dimension)
fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim2);

// reorders any two dimensional integer array to the new mapping (1D mapping must correspond to slowest growing dimension)
int* reorder_2D_int_arr(int* var, int* mapping, int dim1, int dim2);

GSMD_struct reorder_GSMD_struct(int* mapping, GSMD_struct GSMD);

#endif
