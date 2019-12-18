// Written by Samuel Elliott, Summer 2017.

#include <matrix_transformations.h>

#include <stdlib.h>
#include <stdio.h>

// pads 2D fType valued variable and fills with the boundary values
fType* pad_2D_fp_arr(fType* var, int dim1, int padded_dim1, int dim2, int padded_dim2) {

    if (dim1 == padded_dim1 && dim2 == padded_dim2)
        return(var);

    // allocate padded variable
    fType* padded_var = (fType*) malloc(sizeof(fType) * padded_dim1 * padded_dim2);

    // copy var data to padded_var and fill padding with boundary extensions
    for (int i = 0; i < padded_dim1; i++) {
        for (int j = 0; j < padded_dim2; j++) {
            int i_old = i < dim1 ? i : dim1 - 1;	// use boundary row id if in padded row
            int j_old = j < dim2 ? j : dim2 - 1;	// use boundary column id if in padded column
            int k_old = i_old * dim2 + j_old;		// linear index
            int k = i * padded_dim2 + j;			// new linear index
            padded_var[k] = var[k_old];	        	// assign corresponding element
        }
    }

    // free old space and return the padded version
    free(var);

    return padded_var;
}

// pads 2D int valued variable and fills with the boundary values
int* pad_2D_int_arr(int* var, int dim1, int padded_dim1, int dim2, int padded_dim2) {

    if (dim1 == padded_dim1 && dim2 == padded_dim2)
        return(var);

    // allocate padded variable
    int* padded_var = (int*) malloc(sizeof(int) * padded_dim1 * padded_dim2);

    // copy var data to padded_var and fill padding with boundary extensions
    for (int i = 0; i < padded_dim1; i++) {
        for (int j = 0; j < padded_dim2; j++) {
            int i_old = i < dim1 ? i : dim1 - 1;	// use boundary row id if in padded row
            int j_old = j < dim2 ? j : dim2 - 1;	// use boundary column id if in padded column
            int k_old = i_old * dim2 + j_old;		// linear index
            int k = i * padded_dim2 + j;			// new linear index
            padded_var[k] = var[k_old];		// assign corresponding element
        }
    }

    // free old space and return the padded version
    free(var);

    return padded_var;
}

// transposes 2d floating point array
fType* transpose_2D_fp_arr(fType* arr, int dim1, int dim2) {

    fType* transposed_arr = (fType*) malloc(sizeof(fType) * dim1 * dim2);

    for (int i = 0; i < dim2; i++) {
        for (int j = 0; j < dim1; j++) {
            transposed_arr[(i * dim1) + j] = arr[(j * dim2) + i];
        }
    }

    free(arr);

    return transposed_arr;
}

// transposes 2d int array
int* transpose_2D_int_arr(int* arr, int dim1, int dim2) {

    int* transposed_arr = (int*) malloc(sizeof(int) * dim1 * dim2);

    for (int i = 0; i < dim2; i++) {
        for (int j = 0; j < dim1; j++) {
            transposed_arr[(i * dim1) + j] = arr[(j * dim2) + i];
        }
    }

    free(arr);

    return transposed_arr;
}

// tile 2d floating point array
fType* tile_2D_fp_arr(fType* arr, int dim1, int dim2, int dim2_tiled) {
	
    // allocate space for tiled array
    fType* tiled_arr = (fType*) malloc(sizeof(fType) * dim1 * dim2);

    // get number of tiles, tile size, etc.
    int n_tiles = dim2 / dim2_tiled;
    int tile_size = dim1 * dim2_tiled;

    for (int i = 0; i < n_tiles; i++) {		// iterate through each tile
        for (int j = 0; j < dim1; j++) {		// iterate through each row
            for (int k = 0; k < dim2_tiled; k++) {		// iterate through each column
                tiled_arr[(i * tile_size) + (j * dim2_tiled) + k] = arr[ (j * dim2) + (i * dim2_tiled) + k];
            }
        }
    }

    free(arr);

    return tiled_arr;
}

// tile 2d int array
int* tile_2D_int_arr(int* arr, int dim1, int dim2, int dim2_tiled) {
	
    // allocate space for tiled array
    int* tiled_arr = (int*) malloc(sizeof(int) * dim1 * dim2);

    // get number of tiles, tile size, etc.
    int n_tiles = dim2 / dim2_tiled;
    int tile_size = dim1 * dim2_tiled;

    for (int i = 0; i < n_tiles; i++) {		// iterate through each tile
        for (int j = 0; j < dim1; j++) {		// iterate through each row
            for (int k = 0; k < dim2_tiled; k++) {		// iterate through each column
                tiled_arr[(i * tile_size) + (j * dim2_tiled) + k] = arr[ (j * dim2) + (i * dim2_tiled) + k];
            }
        }
    }

    free(arr);

    return tiled_arr;
}

void print_2D_fp_arr(fType* arr, int n1, int n2, int n1_s, int n2_s, int n1_e, int n2_e) {

    printf("\n\n===================================== PRINTING SUBMATRIX =================================\n\n");

    for (int i = n1_s; i < n1_e; i++) {
        printf("\n");
        for (int j = n2_s; j < n2_e; j++) {
            printf("%.8f, ", arr[(i*n2) + j]);
        }
    }
}

void print_2D_int_arr(int* arr, int n1, int n2, int n1_s, int n2_s, int n1_e, int n2_e) {

    printf("\n\n===================================== PRINTING SUBMATRIX =================================\n\n");

    for (int i = n1_s; i < n1_e; i++) {
        printf("\n");
        for (int j = n2_s; j < n2_e; j++) {
            printf("%d, ", arr[(i*n2) + j]);
        }
    }
}
