#ifndef _REORDER_NODES_H_
#define _REORDER_NODES_H_
#include <floating_types.h>

// reorders any one dimensional floating point array according to a one-to-one mapping

fType* reorder_1D_fp_arr(fType* var, int* mapping, int dim);

// reorders any one dimensional integer array according to a one-to-one mapping

int* reorder_1D_int_arr(int* var, int* mapping, int dim);

// reorders any two dimensional floating point array [dim1][dim2] where dim2 is padded to length dim1_stride
// to the new mapping (1D mapping must correspond to dim1, the slowest growing dimension)

fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim1_stride, int dim2_count);

// reorders any two dimensional integer array [dim1][dim2] where dim2 is padded to length dim1_stride
// to the new mapping (1D mapping must correspond to dim1, the slowest growing dimension)

int* reorder_2D_int_arr(int* var, int* mapping, int dim1, int dim1_stride, int dim2_count);

#endif
