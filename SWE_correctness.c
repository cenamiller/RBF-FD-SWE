#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vmath.h>

#include <floating_types.h>
#include <SWE.h>

// Compare S to S after 100 timesteps

void SWE_correctness(SWE_struct SWE, char *inputFile, fType* S, fType* S_ref){

#if 0
  const int U_SV_ID = 0;
  const int V_SV_ID = 1;
  const int W_SV_ID = 2;
  const int H_SV_ID = 3;
#endif

  // Check answers

  fType *du = malloc(sizeof(fType) * SWE.NNodes);
  fType *dv = malloc(sizeof(fType) * SWE.NNodes);
  fType *dw = malloc(sizeof(fType) * SWE.NNodes);
  fType *dh = malloc(sizeof(fType) * SWE.NNodes);

  for(int i=0; i< SWE.NNodes; i++){
    int iu = SWE.NState*i+U_SV_ID;
    int iv = SWE.NState*i+V_SV_ID;
    int iw = SWE.NState*i+W_SV_ID;
    int ih = SWE.NState*i+H_SV_ID;

    du[i] = fabs(S[iu] - S_ref[iu]);
    dv[i] = fabs(S[iv] - S_ref[iv]);
    dw[i] = fabs(S[iw] - S_ref[iw]);
    dh[i] = fabs(S[ih] - S_ref[ih]);
  }

  fType dumax = fmaxval(du, SWE.NNodes);
  fType dvmax = fmaxval(dv, SWE.NNodes);
  fType dwmax = fmaxval(dw, SWE.NNodes);
  fType dhmax = fmaxval(dh, SWE.NNodes);

  printf("\n==============================================================\n\n");
  printf("L1 norm errors MAX(|S_obs - S_ref|) :\n\n");
  printf("u L1 error = %e\n", dumax);
  printf("v L1 error = %e\n", dvmax);
  printf("w L1 error = %e\n", dwmax);
  printf("h L1 error = %e\n", dhmax);
  printf("\n==============================================================\n");

  free(du);
  free(dv);
  free(dw);
  free(dh);
}
