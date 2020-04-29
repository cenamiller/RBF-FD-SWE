#include <stdio.h>
#include <SWE.h>

// Only need _VERBOSE_ to debug initialization
#undef _VERBOSE_

void SWE_init_testcase5(SWE_struct *SWE, fType *S, int *mapping){

  const fType u0 = 20.0;    // Test case 5 zonal velocity 
  //  const fType h0 = 5400.0;
  const int NState = SWE->NState;

  fType *x  = SWE->x;
  fType *y  = SWE->y;
  fType *z  = SWE->z;
  int iu,iv,iw,ih;

  for(int i = 0; i<SWE->NNodes; i++){
    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    S[iu] = -u0*y[i];
    S[iv] = u0*x[i];
    S[iw] = 0.0;
    S[ih] = - (SWE->a*SWE->Omega*u0 + 0.5*u0*u0)*z[i]*z[i];
#ifdef _VERBOSE_
    if (i==0){
      printf("in init tc5: %f %f %f %f\n",S[iu],S[iv],S[iw],S[ih]);
    }
#endif
  }

}

