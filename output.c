#include <stdlib.h>
#include <stdio.h>
#include <floating_types.h>

void output(const fType t, const int NState, const fType *S){
  for(int i=0; i<NState; i++){
    printf("%.15f ",S[i]);
  }
  printf("\n");
}
