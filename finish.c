#include <stdlib.h>
#include <floating_types.h>

void finish(fType *S, fType *S_100){
  free(S);
  free(S_100);
}
