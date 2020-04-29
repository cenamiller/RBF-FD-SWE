#ifndef _SWE_
#define _SWE_
#include <floating_types.h>

#define U_SV_ID 0
#define V_SV_ID 1
#define W_SV_ID 2
#define H_SV_ID 3

typedef struct SWE_struct {

  // --------------- Nodeset and Relevant Nodepoint Data ---------------------------------- //

  int NNodes;        // number of nodes
  int NNbrs;         // number of neighbors
  int NNbrs_Padded;  // total number of nodes padded by some power of two
  int NState;        // Number of state variables in SWEs (e.g. 4)

  // static node coordinates

  fType a;           // radius of Earth (meter)
  fType gh0;         // reference geopotential (meter^2/sec^2)
  fType dt;          // timestep (sec)
  fType gamma;       // diffusion coefficient (sec^-1)
  fType g;           // Gravity acceleration at Earth's surfce (omitted in init file)
  fType Omega;       // angular velocity of Earth (sec^-1)

  fType* x; // x-coordinates
  fType* y; // y-coordinates
  fType* z; // z-coordinates

  fType* f;        // Coriolis force
  fType* ghm;      // the profile of the mountain

  // surface projection operators

  fType* p_u;
  fType* p_v;
  fType* p_w;

  fType* gradghm;  // gradient of the mountains profile

  // ------------------- RBF-FD Differentiation Matrix Data ------------------------------- //

  int* idx;       // Neighbor table
  fType* Dx;      // d/dx weights
  fType* Dy;      // d/dy weights
  fType* Dz;      // d/dz weights
  fType* L;       // Laplacian weights
  
} SWE_struct;

void SWE_reorder_nodes(SWE_struct *SWE,int *mapping);

#endif
