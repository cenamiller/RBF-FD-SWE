#include <stdio.h>
#include <math.h>
#include <floating_types.h>

void fsimple( const fType t, const fType* S, fType* F){

  fType y;

  y=S[0];
  F[0]=(5.0*t*t-y)/exp(t+y);

}

void faero( const fType t, const fType* S, fType* F){

  const fType alpha = 0.098; //  b/m
  const fType g = 9.8;       //  gravity 

  fType x = S[0];
  fType u = S[1];
   
  F[0] = u;
  F[1] = alpha*u*u - g;

}

#include <SWE.h>
extern SWE_struct* SWE;

#include "SWE.h"
#include "Vec3D.h"

void fswe( const fType t, const fType* restrict S, fType* restrict K){

#if 0
  const int X_COMP_ID = 0;
  const int Y_COMP_ID = 1;
  const int Z_COMP_ID = 2;
#endif

  // Dimensions

  const int NNodes = SWE->NNodes;
  const int NState = SWE->NState;
  const int NNbrs  = SWE->NNbrs;
  const int NNbrs_Padded = SWE->NNbrs_Padded;

  // 3-D node coordinates (x,y,z)

  const fType* x = SWE->x;
  const fType* y = SWE->y;
  const fType* z = SWE->z;

  // spherical projetion operators 

  const fType* p_u = SWE->p_u;
  const fType* p_v = SWE->p_v;
  const fType* p_w = SWE->p_w;

  // planetary variables

  const fType  gh0 = SWE->gh0;          // reference geopotential (gh0)
  const fType* f = SWE->f;              // coreolis term
  const fType* ghm = SWE->ghm;          // mountain geopotential height
  const fType* gradghm = SWE->gradghm;  // gradient of geopotential height

  // stencil and derivative operators

  const int*  idx = SWE->idx;
  const fType* Dx = SWE->Dx;
  const fType* Dy = SWE->Dy;
  const fType* Dz = SWE->Dz;
  const fType* L  = SWE->L;

  for(int i=0; i<NNodes; i++){

    fType u;
    fType v;
    fType w;
    fType h;

    // initialize temporary variables for spatial differentiation at respective node

    fType du_dx = 0.0;
    fType dv_dx = 0.0;
    fType dw_dx = 0.0;
    fType dh_dx = 0.0;

    fType du_dy = 0.0;
    fType dv_dy = 0.0;
    fType dw_dy = 0.0;
    fType dh_dy = 0.0;

    fType du_dz = 0.0;
    fType dv_dz = 0.0;
    fType dw_dz = 0.0;
    fType dh_dz = 0.0;
    
    fType Lu = 0.0;
    fType Lv = 0.0;
    fType Lw = 0.0;
    fType Lh = 0.0;

    int iu;
    int iv;
    int iw;
    int ih;

    for (int j = 0; j < NNbrs; j++) {

      int dm_lind = i*NNbrs_Padded + j;

      int inbr = idx[dm_lind];

      iu = U_SV_ID + NState*inbr;
      iv = V_SV_ID + NState*inbr;
      iw = W_SV_ID + NState*inbr;
      ih = H_SV_ID + NState*inbr;

    // temporary variables to hold state variable data for current node

      fType u = S[iu];
      fType v = S[iv];
      fType w = S[iw];
      fType h = S[ih];

      // get neighbor node's state vars

      // update sums

      du_dx += Dx[dm_lind] * u;
      dv_dx += Dx[dm_lind] * v;
      dw_dx += Dx[dm_lind] * w;
      dh_dx += Dx[dm_lind] * h;

      du_dy += Dy[dm_lind] * u;
      dv_dy += Dy[dm_lind] * v;
      dw_dy += Dy[dm_lind] * w;
      dh_dy += Dy[dm_lind] * h;
    
      du_dz += Dz[dm_lind] * u;
      dv_dz += Dz[dm_lind] * v;
      dw_dz += Dz[dm_lind] * w;
      dh_dz += Dz[dm_lind] * h;
        
      Lu += L[dm_lind] * u;
      Lv += L[dm_lind] * v;
      Lw += L[dm_lind] * w;
      Lh += L[dm_lind] * h;
    }

    iu = U_SV_ID + NState*i;
    iv = V_SV_ID + NState*i;
    iw = W_SV_ID + NState*i;
    ih = H_SV_ID + NState*i;

    u = S[iu];
    v = S[iv];
    w = S[iw];
    h = S[ih];

    fType rhs_u = - ((u * du_dx) + (v * du_dy) + (w * du_dz) + (f[i] * ((y[i] * w) - (z[i] * v))) + dh_dx);
    fType rhs_v = - ((u * dv_dx) + (v * dv_dy) + (w * dv_dz) + (f[i] * ((z[i] * u) - (x[i] * w))) + dh_dy);
    fType rhs_w = - ((u * dw_dx) + (v * dw_dy) + (w * dw_dz) + (f[i] * ((x[i] * v) - (y[i] * u))) + dh_dz);

    // vector linear ids

    int vect_idx = 3*i + X_COMP_ID;
    int vect_idy = 3*i + Y_COMP_ID;
    int vect_idz = 3*i + Z_COMP_ID;

    // evaluate projections and apply hyperviscosity

    K[iu] = (p_u[vect_idx] * rhs_u) + (p_u[vect_idy] * rhs_v) + (p_u[vect_idz] * rhs_w) + Lu;
    K[iv] = (p_v[vect_idx] * rhs_u) + (p_v[vect_idy] * rhs_v) + (p_v[vect_idz] * rhs_w) + Lv;
    K[iw] = (p_w[vect_idx] * rhs_u) + (p_w[vect_idy] * rhs_v) + (p_w[vect_idz] * rhs_w) + Lw;

    // ----------------- Use RBF-FD Approximations to Calculate the RHS of the Geopotential Equation ---------------------//

    K[ih] = - ((u * (dh_dx - gradghm[vect_idx])) + (v * (dh_dy - gradghm[vect_idy])) + (w * (dh_dz - gradghm[vect_idz])) 
	    + ((h + gh0 - ghm[i]) * (du_dx + dv_dy + dw_dz))) + Lh;
  }

}


