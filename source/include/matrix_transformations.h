#ifndef MATRIX_TRANSFORMATIONS_H
#define MATRIX_TRANSFORMATIONS_H

#include <swe_config.h>

// pads 2D variable and fills with the boundary values
fType* pad_2D_fp_arr(fType* var, int dim1, int padded_dim1, int dim2, int padded_dim2);

// pads 2D int valued variable and fills with the boundary values
int* pad_2D_int_arr(int* var, int dim1, int padded_dim1, int dim2, int padded_dim2);

fType* transpose_2D_fp_arr(fType* arr, int dim1, int dim2);

int* transpose_2D_int_arr(int* arr, int dim1, int dim2);

// tile 2d floating point array
fType* tile_2D_fp_arr(fType* arr, int dim1, int dim2, int dim2_tiled);	

// tile 2d int array
int* tile_2D_int_arr(int* arr, int dim1, int dim2, int dim2_tiled);

void print_2D_fp_arr(fType* arr, int n1, int n2, int n1_s, int n2_s, int n1_e, int n2_e);
void print_2D_int_arr(int* arr, int n1, int n2, int n1_s, int n2_s, int n1_e, int n2_e);

#endif
