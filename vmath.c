#include <floating_types.h>
#include <stdio.h>

fType fminval(fType *a, int vlen){
  fType vmin;
  vmin = a[0];
  for(int i=1; i< vlen; i++){
    if (a[i] < vmin){
      vmin = a[i];
    }
  }   
  return vmin;
}

fType fmaxval(fType *a, int vlen){
  fType vmax;
  vmax = a[0];
  for(int i=1; i< vlen; i++){
    if (a[i] > vmax){
      vmax = a[i];
    }
  }   
  return vmax;
}

fType fsum(fType *a, int vlen){
  fType vsum; 
  vsum = 0.0;
  for(int i=0; i< vlen; i++){
    vsum +=a[i];
  }
  return vsum;
}
