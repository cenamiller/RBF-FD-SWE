#include <stdio.h>
#include <stdlib.h>
#include <floating_types.h>

#define CHECK_OPEN_ERR(f, i) {if (f) { printf("ERROR in SWE_final_conditions: Could not open %s while trying to read data \n", i); exit(2);}}

void SWE_read_initial_conditions(char *inputFile, fType *S){

  long NNodes, NNbrs;
  int temp;

  // RDL Enables fType and fType_Input to actually differ

  double *tmp;

  // open input file
  FILE* fp_in = fopen(inputFile, "r");
  CHECK_OPEN_ERR(!fp_in, inputFile);

  // constants/parameters
  fread((void*) &temp, sizeof(int), 1, fp_in);
  NNodes = (long) temp;
  fread((void*) &temp, sizeof(int), 1, fp_in);
  NNbrs = (long) temp;

  // calculate and position file pointer
  long offset = (sizeof(int) * (2 + (NNodes * NNbrs))) + (sizeof(double) * (4 + (5 * NNodes) + (4 * NNodes * 3) + (4 * NNodes * NNbrs)));
  fseek(fp_in, offset, SEEK_SET);

  // allocate/read H_init
  // RDL Enables fType and fType_Input to actually differ

  tmp = (double*) malloc(sizeof(double) * NNodes * 4);
  fread((void*) tmp, sizeof(double), NNodes * 4, fp_in);
  for(int i=0; i<4*NNodes; i++){
    S[i] = (fType) tmp[i];
  }
  
  free(tmp);

  fclose(fp_in);
}
