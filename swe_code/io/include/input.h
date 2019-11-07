#ifndef SWE_INPUT_H
#define SWE_INPUT_H

#include <swe_config.h>

#ifdef USE_NCIO

// read simulation GSMD struct data
GSMD_struct get_GSMD_nc(char* inputFile);

// read initial conditions for the state variable matrix
fType* get_ICs_nc(char* inputFile);

fType* get_FCs_nc(char* inputFile);

#endif

GSMD_struct get_GSMD_bin(char* inputFile);

fType* get_ICs_bin(char* inputFile);

fType* get_FCs_bin(char* inputFile);

#endif
