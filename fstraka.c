#include <NS.h>
extern NS_struct* NS;

void fstraka( const fType t, const fType* restrict S, fType* restrict K){

  const int U_SV_ID  = 0;
  const int W_SV_ID  = 1;
  const int TH_SV_ID = 2;
  const int PI_SV_ID = 3;

  // Dimensions

  const int NNodes = NS->NNodes;
  const int NState = NS->NState;
  const int NNbrs  = NS->NNbrs;
  const int NNbrs_Padded = NS->NNbrs_Padded;

  // planetary variables

  const fType  g = NS->g;          //
  const fType  Cp = NS->Cp;        //
  const fType  Cv = NS->Cv;        // 
  const fType  Rd = NS->Rd;        //
  const fType  mu = NS->mu;        //
  const fType  thbar = NS->Ts;      // background temperature from Straka test case

  fType dpibar_dz = - g/(Cp*thbar);
  fType *pibar = NS->pibar;

  // stencil and derivative operators

  const int*  idx = NS->idx;
  const fType* Dx = NS->Dx;
  const fType* Dz = NS->Dz;
  const fType* L  = NS->L;

  for(int i=0; i<NNodes; i++){

    // initialize temporary variables for spatial differentiation at respective node

    fType u  = 0.0;
    fType w  = 0.0;
    fType th  = 0.0;
    fType pi  = 0.0;

    fType du_dx  = 0.0;
    fType dw_dx  = 0.0;
    fType dth_dx = 0.0;
    fType dpi_dx = 0.0;

    fType du_dz  = 0.0;
    fType dw_dz  = 0.0;
    fType dth_dz = 0.0;
    fType dpi_dz = 0.0;
    
    fType Lu = 0.0;
    fType Lw = 0.0;
    fType Lth = 0.0;

    for (int j = 0; j < NNbrs; j++) {

      int dm_lind = i*NNbrs_Padded + j;

      int inbr = idx[dm_lind];

      int iu  = U_SV_ID + NState*inbr;
      int iw  = W_SV_ID + NState*inbr;
      int ith = TH_SV_ID + NState*inbr;
      int ipi = PI_SV_ID + NState*inbr;

      u  = S[iu];
      w  = S[iw];
      th = S[ith];
      pi = S[ipi];

      // get neighbor node's state vars

      // update sums

      du_dx  += Dx[dm_lind] * u;
      dw_dx  += Dx[dm_lind] * w;
      dpi_dx += Dx[dm_lind] * pi;
      dth_dx += Dx[dm_lind] * th;

      du_dz  += Dz[dm_lind] * u;
      dw_dz  += Dz[dm_lind] * w;
      dpi_dz += Dz[dm_lind] * pi;
      dth_dz += Dz[dm_lind] * th;
        
      Lu  += L[dm_lind] * u;
      Lw  += L[dm_lind] * w;
      Lth += L[dm_lind] * th;
    }

    int iu  = U_SV_ID + NState*i;
    int iw  = W_SV_ID + NState*i;
    int ith = TH_SV_ID + NState*i;
    int ipi = PI_SV_ID + NState*i;

    u  = S[iu];
    w  = S[iw];
    th = S[ith];
    pi = S[ipi];

    fType rhs_u  = - u * du_dx  - w * du_dz - Cp*(thbar+th)*dpi_dx + mu*Lu;
    fType rhs_w  = - u * dw_dx  - w * du_dz - Cp*(thbar+th)*dpi_dz + (g/thbar)*th + mu*Lw;
    fType rhs_th = - u * dth_dx - w * dth_dz +  mu*Lth;
    fType rhs_pi = - u * dpi_dx - w * (dpibar_dz + dpi_dz) - (Rd/Cv)*(pibar[i]+pi)*(du_dx + dw_dz);

    // ----------------------------------------------------

    // evaluate projections and apply hyperviscosity

    K[iu]  = rhs_u;
    K[iw]  = rhs_w;
    K[ith] = rhs_th;
    K[ipi] = rhs_pi;
  }

}

