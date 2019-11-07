// Written by Samuel Elliot, Summer 2017.

#include <reorder_nodes.h>

#include <stdlib.h>

// reorders any two dimensional floating point array to the new mapping (1D mapping must correspond to slowest growing dimension)
fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim2) {
	
    // create space for mapped variable
    fType* mapped_var = (fType *) malloc(sizeof(fType) * dim1 * dim2);

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            mapped_var[(i * dim2) + j] = var[(mapping[i] * dim2) + j];
        }
    }

    // free memory allocated to the initial variable data
    free(var);

    // return the new variable data
    return mapped_var;
}

// reorders any two dimensional integer array to the new mapping (1D mapping must correspond to slowest growing dimension)
int* reorder_2D_int_arr(int* var, int* mapping, int dim1, int dim2) {
	
    // create space for mapped variable
    int* mapped_var = (int *) malloc(sizeof(int) * dim1 * dim2);

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            mapped_var[(i * dim2) + j] = var[(mapping[i] * dim2) + j];
        }
    }

    // free memory allocated to the initial variable data
    free(var);

    // return the new variable data
    return mapped_var;
}


GSMD_struct reorder_GSMD_struct(int* mapping, GSMD_struct GSMD) {

    // replace all variable data pointed to by the GSMD struct with reordered versions
    // scalar variables
    GSMD.x = reorder_2D_fp_arr(GSMD.x, mapping, GSMD.Nnodes, 1);
    GSMD.y = reorder_2D_fp_arr(GSMD.y, mapping, GSMD.Nnodes, 1);
    GSMD.z = reorder_2D_fp_arr(GSMD.z, mapping, GSMD.Nnodes, 1);
    GSMD.f = reorder_2D_fp_arr(GSMD.f, mapping, GSMD.Nnodes, 1);
    GSMD.ghm = reorder_2D_fp_arr(GSMD.ghm, mapping, GSMD.Nnodes, 1);

    // vector variables
    GSMD.p_u = reorder_2D_fp_arr(GSMD.p_u, mapping, GSMD.Nnodes, 3);
    GSMD.p_v = reorder_2D_fp_arr(GSMD.p_v, mapping, GSMD.Nnodes, 3);
    GSMD.p_w = reorder_2D_fp_arr(GSMD.p_w, mapping, GSMD.Nnodes, 3);
    GSMD.gradghm = reorder_2D_fp_arr(GSMD.gradghm, mapping, GSMD.Nnodes, 3);

    // reorder all weights
    GSMD.Dx = reorder_2D_fp_arr(GSMD.Dx, mapping, GSMD.Nnodes, GSMD.Nnbr);
    GSMD.Dy = reorder_2D_fp_arr(GSMD.Dy, mapping, GSMD.Nnodes, GSMD.Nnbr);
    GSMD.Dz = reorder_2D_fp_arr(GSMD.Dz, mapping, GSMD.Nnodes, GSMD.Nnbr);
    GSMD.L = reorder_2D_fp_arr(GSMD.L, mapping, GSMD.Nnodes, GSMD.Nnbr);

    // reorder idx
    GSMD.idx = reorder_2D_int_arr(GSMD.idx, mapping, GSMD.Nnodes, GSMD.Nnbr);

    // create inverse mapping: inv_mapping[old_idx] == new_idx
    int* inv_mapping = (int*) malloc(sizeof(int)*GSMD.Nnodes);

    for (int i = 0; i < GSMD.Nnodes; i++)
        inv_mapping[mapping[i]] = i;

    // account for node id reordering in idx
    for (int i = 0; i < GSMD.Nnodes * GSMD.Nnbr; i++) {
        GSMD.idx[i] = inv_mapping[GSMD.idx[i]];
    }

    // free inverse mapping data
    free(inv_mapping);

    // return reordered GSMD
    return GSMD;
}
