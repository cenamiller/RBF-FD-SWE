// Written by Samuel Elliott, Summer 2017 (code is unused in project but helps convert nc to bin filetypes)

#include <swe_config.h>
#include <input.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
	
    // file configuration
    char inputFile[MAX_FILE_PATH_SIZE];
    char outputFile[MAX_FILE_PATH_SIZE];
    strcpy(&inputFile[0],getenv("NC2BIN_INPUT_FILE"));
    strcpy(&outputFile[0],getenv("NC2BIN_OUTPUT_FILE"));

    // Data
    fType* H_init;
    fType* H_100;
    GSMD_struct GSMD;

    // read inputs
    GSMD = get_GSMD_nc(inputFile);
    H_init = get_ICs_nc(inputFile);
    H_100 = get_FCs_nc(inputFile);

    // open output file
    FILE* fp_out = fopen(outputFile, "w+");

    // constants/parameters
    fwrite((void*) &GSMD.Nnodes, sizeof(int), 1, fp_out);
    fwrite((void*) &GSMD.Nnbr, sizeof(int), 1, fp_out);
    fwrite((void*) &GSMD.a, sizeof(double), 1, fp_out);
    fwrite((void*) &GSMD.gh0, sizeof(double), 1, fp_out);
    fwrite((void*) &GSMD.dt, sizeof(double), 1, fp_out);
    fwrite((void*) &GSMD.gamma, sizeof(double), 1, fp_out);

    printf("\nNC2BIN:\tNnodes = %d,\tNnbr = %d\n",GSMD.Nnodes,GSMD.Nnbr);
    
    // scalar vars
    fwrite((void*) GSMD.x, sizeof(double), GSMD.Nnodes, fp_out);
    fwrite((void*) GSMD.y, sizeof(double), GSMD.Nnodes, fp_out);
    fwrite((void*) GSMD.z, sizeof(double), GSMD.Nnodes, fp_out);
    fwrite((void*) GSMD.f, sizeof(double), GSMD.Nnodes, fp_out);
    fwrite((void*) GSMD.ghm, sizeof(double), GSMD.Nnodes, fp_out);

    // vector vars
    fwrite((void*) GSMD.p_u, sizeof(double), GSMD.Nnodes * 3, fp_out);
    fwrite((void*) GSMD.p_v, sizeof(double), GSMD.Nnodes * 3, fp_out);
    fwrite((void*) GSMD.p_w, sizeof(double), GSMD.Nnodes * 3, fp_out);
    fwrite((void*) GSMD.gradghm, sizeof(double), GSMD.Nnodes * 3, fp_out);

    // unscale differentiation weights
    for (int i = 0; i < GSMD.Nnodes * GSMD.Nnbr; i++) {
        GSMD.Dx[i] = GSMD.Dx[i]*GSMD.a;
        GSMD.Dy[i] = GSMD.Dy[i]*GSMD.a;
        GSMD.Dz[i] = GSMD.Dz[i]*GSMD.a;
        GSMD.L[i] = GSMD.L[i]/GSMD.gamma;
    }

    // write DMS
    fwrite((void*) GSMD.Dx, sizeof(double), GSMD.Nnodes * GSMD.Nnbr, fp_out);
    fwrite((void*) GSMD.Dy, sizeof(double), GSMD.Nnodes * GSMD.Nnbr, fp_out);
    fwrite((void*) GSMD.Dz, sizeof(double), GSMD.Nnodes * GSMD.Nnbr, fp_out);
    fwrite((void*) GSMD.L, sizeof(double), GSMD.Nnodes * GSMD.Nnbr, fp_out);

    // unshift stencil indices
    for (int i = 0; i < GSMD.Nnodes * GSMD.Nnbr; i++) {
        GSMD.idx[i] = GSMD.idx[i] + 1;
    }

    // stencils
    fwrite((void*) GSMD.idx, sizeof(int), GSMD.Nnodes * GSMD.Nnbr, fp_out);

    // initial conditions
    fwrite((void*) H_init, sizeof(double), GSMD.Nnodes * 4, fp_out);

    // solution at 100 timesteps
    fwrite((void*) H_100, sizeof(double), GSMD.Nnodes * 4, fp_out);

    // close binary output file
    fclose(fp_out);
}
