// Written by Samuel Elliott, Summer 2017.

#include <layout.h>
#include <matrix_transformations.h>

extern sim_params_struct sim_params;

    // pads all GSMD data to padded_Nnodes which aligns with PADDED_NUM
    // fills new nodes with the same data as the last node -> simply recalculates these values and doesn't cause any parallelization/distribution issues
    GSMD_struct pad_GSMD_data(GSMD_struct GSMD) {

    // get padded_nnodes
    GSMD.padded_Nnodes = ROUNDUP_SIMD(GSMD.Nnodes);

    // pad all the corresponding GSMD data
    // coordinates
    GSMD.x = pad_2D_fp_arr(GSMD.x, GSMD.Nnodes, GSMD.padded_Nnodes, 1, 1);
    GSMD.y = pad_2D_fp_arr(GSMD.y, GSMD.Nnodes, GSMD.padded_Nnodes, 1, 1);
    GSMD.z = pad_2D_fp_arr(GSMD.z, GSMD.Nnodes, GSMD.padded_Nnodes, 1, 1);

    // scalars
    GSMD.f = pad_2D_fp_arr(GSMD.f, GSMD.Nnodes, GSMD.padded_Nnodes, 1, 1);
    GSMD.ghm = pad_2D_fp_arr(GSMD.ghm, GSMD.Nnodes, GSMD.padded_Nnodes, 1, 1);

    // spacial vectors
    GSMD.p_u = pad_2D_fp_arr(GSMD.p_u, GSMD.Nnodes, GSMD.padded_Nnodes, 3, 3);
    GSMD.p_v = pad_2D_fp_arr(GSMD.p_v, GSMD.Nnodes, GSMD.padded_Nnodes, 3, 3);
    GSMD.p_w = pad_2D_fp_arr(GSMD.p_w, GSMD.Nnodes, GSMD.padded_Nnodes, 3, 3);
    GSMD.gradghm = pad_2D_fp_arr(GSMD.gradghm, GSMD.Nnodes, GSMD.padded_Nnodes, 3, 3);

    #ifdef USE_SFDL
    GSMD.padded_Nnbr = GSMD.Nnbr;
    #else
    GSMD.padded_Nnbr = ROUNDUP_SIMD(GSMD.Nnbr);
    #endif

    // pad Dx,Dy,Dz,L and idx
    GSMD.Dx = pad_2D_fp_arr(GSMD.Dx, GSMD.Nnodes, GSMD.padded_Nnodes, GSMD.Nnbr, GSMD.padded_Nnbr);
    GSMD.Dy = pad_2D_fp_arr(GSMD.Dy, GSMD.Nnodes, GSMD.padded_Nnodes, GSMD.Nnbr, GSMD.padded_Nnbr);
    GSMD.Dz = pad_2D_fp_arr(GSMD.Dz, GSMD.Nnodes, GSMD.padded_Nnodes, GSMD.Nnbr, GSMD.padded_Nnbr);
    GSMD.L = pad_2D_fp_arr(GSMD.L, GSMD.Nnodes, GSMD.padded_Nnodes, GSMD.Nnbr, GSMD.padded_Nnbr);
    GSMD.idx = pad_2D_int_arr(GSMD.idx, GSMD.Nnodes, GSMD.padded_Nnodes, GSMD.Nnbr, GSMD.padded_Nnbr);

    return GSMD;
}


GSMD_struct convert2_SFDL(GSMD_struct GSMD) {

    // transpose 3D vector data
    GSMD.p_u = transpose_2D_fp_arr(GSMD.p_u, GSMD.padded_Nnodes, 3);
    GSMD.p_v = transpose_2D_fp_arr(GSMD.p_v, GSMD.padded_Nnodes, 3);
    GSMD.p_w = transpose_2D_fp_arr(GSMD.p_w, GSMD.padded_Nnodes, 3);
    GSMD.gradghm = transpose_2D_fp_arr(GSMD.gradghm, GSMD.padded_Nnodes, 3);

    // tile 3D vector data
    GSMD.p_u = tile_2D_fp_arr(GSMD.p_u, 3, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.p_v = tile_2D_fp_arr(GSMD.p_v, 3, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.p_w = tile_2D_fp_arr(GSMD.p_w, 3, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.gradghm = tile_2D_fp_arr(GSMD.gradghm, 3, GSMD.padded_Nnodes, SIMD_LENGTH);

    // first transpose all DMs to get NNMI layout
    GSMD.Dx = transpose_2D_fp_arr(GSMD.Dx, GSMD.padded_Nnodes, GSMD.padded_Nnbr);
    GSMD.Dy = transpose_2D_fp_arr(GSMD.Dy, GSMD.padded_Nnodes, GSMD.padded_Nnbr);
    GSMD.Dz = transpose_2D_fp_arr(GSMD.Dz, GSMD.padded_Nnodes, GSMD.padded_Nnbr);
    GSMD.L = transpose_2D_fp_arr(GSMD.L, GSMD.padded_Nnodes, GSMD.padded_Nnbr);
    GSMD.idx = transpose_2D_int_arr(GSMD.idx, GSMD.padded_Nnodes, GSMD.padded_Nnbr);

    // tile the GSMD with tiles of dimension Nnbr by SIMD_LENGTH
    GSMD.Dx = tile_2D_fp_arr(GSMD.Dx, GSMD.padded_Nnbr, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.Dy = tile_2D_fp_arr(GSMD.Dy, GSMD.padded_Nnbr, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.Dz = tile_2D_fp_arr(GSMD.Dz, GSMD.padded_Nnbr, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.L = tile_2D_fp_arr(GSMD.L, GSMD.padded_Nnbr, GSMD.padded_Nnodes, SIMD_LENGTH);
    GSMD.idx = tile_2D_int_arr(GSMD.idx, GSMD.padded_Nnbr, GSMD.padded_Nnodes, SIMD_LENGTH);

    // return new GSMD struct
    return GSMD;
}
