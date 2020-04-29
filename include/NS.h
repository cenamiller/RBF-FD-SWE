#ifndef _NS_H_
#define _NS_H_
#include <floating_types.h>


typedef struct NS_struct {

  // --------------- Nodeset and Relevant Nodepoint Data ---------------------------------- //

  int NNodes;        // number of nodes
  int NNbrs;         // number of neighbors
  fType dt;          // timestep (sec)

  int NNbrs_Padded;  // total number of nodes padded by some power of two
  int NState;        // Number of state variables in SWEs (e.g. 4)

  // static node coordinates

  fType g;           // Gravity acceleration at Earth's surfce (omitted in init file)
  fType Cp;          // specific heat of air at constant pressure (J/kg)/K
  fType Cv;          // specific heat of air at constant volume (J/kg)/K
  fType Rd;          // gas constant of air (J/kg)/K
  fType mu;          // diffusion coefficient for hyper viscosity (m^2/s)
  fType Ts;          // reference temperature (K)
  fType thbar;       // reference potential temperature (K)
  fType dpibar_dz;   // vertical derivative of Exner pressure

  // vector quantities

  fType* x;          // z coordinate
  fType* z;          // z coordinate
  fType* pibar;      // background Exner pressure

  // ------------------- RBF-FD Differentiation Matrix Data ------------------------------- //

  int*   idx;     // Neighbor table
  fType* Dx;      // d/dx weights
  fType* Dz;      // d/dz weights
  fType* L;       // Laplacian weights
  
} NS_struct;

void NS_reorder_nodes(NS_struct *NS, int *mapping);

#endif
