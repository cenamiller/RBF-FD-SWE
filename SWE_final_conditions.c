#include <stdlib.h>
#include <stdio.h>
#include <floating_types.h>
#include <reorder_nodes.h>

// Get Final Conditions after 100 timesteps

#define CHECK_OPEN_ERR(f, i) {if (f) { printf("ERROR in SWE_final_conditions: Could not open %s while trying to read data \n", i); exit(2);}}

void SWE_final_conditions(char* inputFile, fType* S_ref) {

  long NNodes, NNbrs;
  int temp;

  // RDL Enables fType and fType_Input to actually differ

  double* tmp;

  // open input file

  printf("\n==============================================================\n");
  printf("\nInitializing verification state vector...\n\n");
  printf("Opening file %s",inputFile);

  FILE* fp_in = fopen(inputFile, "r");
  CHECK_OPEN_ERR(!fp_in, inputFile);
  printf(" ...done.\n");

  // constants/parameters
  fread((void*) &temp, sizeof(int), 1, fp_in);
  NNodes = (long) temp;
  fread((void*) &temp, sizeof(int), 1, fp_in);
  NNbrs = (long) temp;
  int NPts = NNodes*4;

  // calculate and position file pointer
  long offset = (sizeof(int) * (2 + (NNodes * NNbrs))) + (sizeof(double) * (4 + (5 * NNodes) + (4 * NNodes * 3) + (4 * NNodes * NNbrs) + (NNodes * 4)));
  fseek(fp_in, offset, SEEK_SET);

  // allocate/read S_ref
    
  tmp = (double*) malloc(sizeof(double) * NPts);
  fread((void*) tmp, sizeof(double), NPts, fp_in);

  for(int i=0; i<NPts; i++){
      S_ref[i] = (fType) tmp[i];
    }

  free(tmp);

  fclose(fp_in);
}
