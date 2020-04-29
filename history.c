#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <floating_types.h>

void history(const char* pdecase, fType *S, fType *S_obs, int NPts){
  if(! strcmp(pdecase, "swe")){
    for(int i=0; i<NPts; i++){
      S_obs[i]=S[i];
    }
  }
  else{
    printf("correctness check not implemented (yet)\n");
  }
}
