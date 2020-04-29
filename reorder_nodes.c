// Written by Samuel Elliot, Summer 2017.

#include <reorder_nodes.h>

#include <stdlib.h>

// reorders any 1 dimensional floating point array to the new mapping
fType* reorder_1D_fp_arr(fType* var, int* mapping, int dim1) {
  
  // create space for mapped variable
  fType* mapped_var = (fType *) malloc(sizeof(fType) * dim1);

  for (int i = 0; i < dim1; i++) {
    mapped_var[i] = var[mapping[i]];
  }

  // free memory allocated to the initial variable data
  free(var);

  // return the new variable data
  return mapped_var;
}

// reorders any 1 dimensional integer array to the new mapping
int* reorder_1D_int_arr(int* var, int* mapping, int dim1) {
  
  // create space for mapped variable
  int* mapped_var = (int *) malloc(sizeof(fType) * dim1);

  for (int i = 0; i < dim1; i++) {
    mapped_var[i] = var[mapping[i]];
  }

  // free memory allocated to the initial variable data
  free(var);

  // return the new variable data
  return mapped_var;
}

// reorders any two dimensional integer array to the new mapping (1D mapping must correspond to slowest growing dimension)
fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim1_stride, int dim2) {
  
  // create space for mapped variable
  fType* mapped_var = (fType *) malloc(sizeof(fType) * dim1 * dim1_stride);

  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      mapped_var[(i * dim1_stride) + j] = var[(mapping[i] * dim1_stride) + j];
    }
  }

  // free memory allocated to the initial variable data
  free(var);

  // return the new variable data
  return mapped_var;

}


// reorders any two dimensional integer array to the new mapping (1D mapping must correspond to slowest growing dimension)
int* reorder_2D_int_arr(int* var, int* mapping, int dim1, int dim1_stride, int dim2) {
  
  // create space for mapped variable
  int* mapped_var = (int *) malloc(sizeof(int) * dim1 * dim1_stride);

  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      mapped_var[(i * dim1_stride) + j] = var[(mapping[i] * dim1_stride) + j];
    }
  }

  // free memory allocated to the initial variable data
  free(var);

  // return the new variable data   
  return mapped_var;
}

