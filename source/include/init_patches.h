#ifndef SWE_INIT_PATCHES_H
#define SWE_INIT_PATCHES_H

#include <swe_config.h>

// determines/assigns global patch decomposition (GPD) for MPI
patch_struct* get_GPD(GSMD_struct GSMD);

void init_LPSMD(patch_struct* GPD, GSMD_struct GSMD);

void verify_valid_GPD(patch_struct* GPD);

fType* get_patch_ICs(fType* H_global, patch_struct* GPD);

void print_patches(patch_struct* GPD);

#endif
