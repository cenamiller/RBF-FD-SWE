#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <floating_types.h>

#include "SWE.h"
#include "Vec3D.h"
#include "vmath.h"

// Only need _VERBOSE_ to debug initialization
#undef _VERBOSE_
#define _TINY_VAL_ 1.e-11

void SWE_init_testcase5(SWE_struct *SWE, fType *S);
void SWE_read_initial_conditions(char *inputFile, fType *S);

#define CHECK_INPUT_ERR(f, in1, in2) {if (f) { printf("ERROR in SWE_init: model value (%i) and file values (%i) do not match!\n", in1, in2); exit(2);}}
#define CHECK_OPEN_ERR(f, i) {if (f) { printf("ERROR in SWE_init: Could not open %s while trying to read data \n", i); exit(2);}}

#define COMPUTE_XYZ  0
#define COMPUTE_FCOR 1
#define COMPUTE_PROJ 1
#define COMPUTE_DERIV 0
#define COMPUTE_TC5 1

// Compare vector from file (vfile) with computed version (vcomp)

fType check_v(double *vfile, fType *vcomp, int NNodes){

  fType *dv = malloc( NNodes * sizeof(fType) );
  for(int i = 0; i < NNodes; i++){
    dv[i] = fabs((fType) vfile[i] - vcomp[i]);
  }
  fType delv = fmaxval(dv, NNodes);
  free(dv);
  return delv;

}

// Compare 3-vector (v3) from file (v3file) with computed version (v3comp)

fType check_v3(double *v3file, fType *v3comp, int NNodes){

  fType *dv3 = malloc( 3*NNodes * sizeof(fType) );
  for(int i = 0; i < 3*NNodes; i++){
    dv3[i] = fabs((fType) v3file[i] - v3comp[i]);
  }
  fType delv3 = fmaxval(dv3, 3*NNodes);
  free(dv3);
  return delv3;
}

// Compare state-vector from file (Sfile) with computed version (Scomp)

fType check_state(fType *Sfile, fType *Scomp, int NState, int NNodes){

  fType *dS = malloc( NState*NNodes * sizeof(fType) );
  for(int i = 0; i < NState*NNodes; i++){
    dS[i] = fabs(Sfile[i] - Scomp[i]);
  }
  fType delS = fmaxval(dS, NState*NNodes);
  free(dS);
  return delS;
}

void SWE_init_bin(char* inputFile, char* testcase, SWE_struct *SWE, fType *S) {

  // RDL temporary variables to enable Ftype and Ftype_Input to differ

  double a;
  double gh0;
  double dt;
  double gamma;

  printf("\n==============================================================\n");
  printf("\nSetting up SWE tc5 initial conditions...\n\n");
  printf("Opening file %s",inputFile);

#if 0
  const int X_COMP_ID = 0;
  const int Y_COMP_ID = 1;
  const int Z_COMP_ID = 2;
#endif

  // open input file

  FILE* fp_in = fopen(inputFile, "r");
  CHECK_OPEN_ERR(!fp_in, inputFile);
  printf(" ...done.\n");

  // constants/parameters

  fread((void*) &SWE->NNodes, sizeof(int), 1, fp_in);
  fread((void*) &SWE->NNbrs, sizeof(int), 1, fp_in);

  //  int NPAD = 8;
  int NPAD = 8;

  SWE->NNbrs_Padded = NPAD*(((SWE->NNbrs-1)/NPAD)+1);
  SWE->NState=4;
  printf("\n\nSWE RBF FD Test configuration:\n\n");
  printf("Stencil size (NNbrs) = %i\nPadded Stencil storage (words) = %i\nNumber of Nodes (NNodes) = %i\n\n",SWE->NNbrs,SWE->NNbrs_Padded,SWE->NNodes);

  fread((void*) &a, sizeof(double), 1, fp_in);
  SWE->a = (fType)a;

  double Omega = 7.292E-5; // Rotational rate of Earth (omitted in init file)
  double g = 9.807;        // Gravity acceleration at Earth's surfce (omitted in init file)

  SWE->Omega = (fType)Omega;
  SWE->g     = (fType)g;

  fread((void*) &gh0, sizeof(double), 1, fp_in);
  SWE->gh0 = (fType)gh0;

  fread((void*) &dt, sizeof(double), 1, fp_in);
  SWE->dt = (fType)dt;

  fread((void*) &gamma, sizeof(double), 1, fp_in);
  SWE->gamma = (fType)gamma;


  printf("\nPhysical parameters:\n\n");
  printf("a      radius of earth (m)           = %e\n",SWE->a);
  printf("g      acceleration of gravity on earth's surface (m/s^2) = %e\n",SWE->g);
  printf("Omega  angular velocity of Earth (sec^-1) = %e\n",SWE->Omega);
  printf("gh0    geopotential height (m^2/s^2) = %e\n",SWE->gh0);
  printf("dt     time step  (s)                = %e\n",SWE->dt);
  printf("gamma  diffusion coeff (s^-1)        = %e\n\n",SWE->gamma);


  // allocate and read coordinates

  double *x = (double*) malloc(sizeof(double) * SWE->NNodes);
  double *y = (double*) malloc(sizeof(double) * SWE->NNodes);
  double *z = (double*) malloc(sizeof(double) * SWE->NNodes);

  SWE->x = (fType*) malloc(sizeof(fType) * SWE->NNodes);
  SWE->y = (fType*) malloc(sizeof(fType) * SWE->NNodes);
  SWE->z = (fType*) malloc(sizeof(fType) * SWE->NNodes);

  fread((void*) x, sizeof(double), SWE->NNodes, fp_in);
  fread((void*) y, sizeof(double), SWE->NNodes, fp_in);
  fread((void*) z, sizeof(double), SWE->NNodes, fp_in);

  if (COMPUTE_XYZ == 1){
    printf("ERROR in SWE_init_bin: Computing (x,y,z) on the sphere not supported\n");
    exit(-1);
  }
  for(int i=0; i<SWE->NNodes; i++){
    SWE->x[i] = (fType) x[i];
    SWE->y[i] = (fType) y[i];
    SWE->z[i] = (fType) z[i];
    }


  // read in coreolis force, or compute it from formula

  double* tmp = (double *) malloc(sizeof(double) * SWE->NNodes);
  fread((void*) tmp, sizeof(double), SWE->NNodes, fp_in);
  SWE->f = (fType*) malloc(sizeof(fType) * SWE->NNodes);
  if (COMPUTE_FCOR == 1){

    // compute f 

    fType* f = (double *) malloc(sizeof(double) * SWE->NNodes);
    for(int i=0; i<SWE->NNodes; i++){
      f[i]   =  2.0*Omega*z[i]; 
      SWE->f[i] = (fType) f[i];
    }
    printf("Computing Coriolis term (f)       ");
    fType f_del = check_v(tmp, SWE->f, SWE->NNodes);
    if (f_del < _TINY_VAL_){
      printf("...PASS\n");
    }
    else{
      printf("...FAIL\n");
      exit(-1);
    }
    free(f);
  }
  else{

    // use the f read from the file

    for(int i=0; i<SWE->NNodes; i++)
      SWE->f[i] = (fType) tmp[i];
  }
  free(tmp);

  double* ghm = (double *) malloc(sizeof(double) * SWE->NNodes);
  SWE->ghm = (fType*) malloc(sizeof(fType) * SWE->NNodes);
  fread((void*) ghm, sizeof(double), SWE->NNodes, fp_in);
  for(int i=0; i<SWE->NNodes; i++)
    SWE->ghm[i] = (fType) ghm[i];
  free(ghm);

  // 3-vector section...

  double* p_u = (double*) malloc(sizeof(double) * SWE->NNodes * 3);
  double* p_v = (double*) malloc(sizeof(double) * SWE->NNodes * 3);
  double* p_w = (double*) malloc(sizeof(double) * SWE->NNodes * 3);

  SWE->p_u  = (fType*) malloc(sizeof(fType) * SWE->NNodes * 3);
  SWE->p_v  = (fType*) malloc(sizeof(fType) * SWE->NNodes * 3);
  SWE->p_w  = (fType*) malloc(sizeof(fType) * SWE->NNodes * 3);

  fread((void*) p_u, sizeof(double), SWE->NNodes * 3, fp_in);
  fread((void*) p_v, sizeof(double), SWE->NNodes * 3, fp_in);
  fread((void*) p_w, sizeof(double), SWE->NNodes * 3, fp_in);

  if (COMPUTE_PROJ == 1){
    for(int i=0; i<SWE->NNodes; i++){
      int idx = i*3 + X_COMP_ID;
      int idy = i*3 + Y_COMP_ID;
      int idz = i*3 + Z_COMP_ID;

      SWE->p_u[idx] = (fType) (1.0 - x[i]*x[i]);
      SWE->p_u[idy] = (fType) (-x[i]*y[i]);
      SWE->p_u[idz] = (fType) (-x[i]*z[i]);

      SWE->p_v[idx] = (fType) (-x[i]*y[i]);
      SWE->p_v[idy] = (fType) (1.0 - y[i]*y[i]);
      SWE->p_v[idz] = (fType) (-y[i]*z[i]);

      SWE->p_w[idx] = (fType) (-x[i]*z[i]);
      SWE->p_w[idy] = (fType) (-y[i]*z[i]);
      SWE->p_w[idz] = (fType) (1.0 - z[i]*z[i]);

    }

    if (COMPUTE_PROJ == 1){
      printf("Computing Projection Operator (u) ");
      fType p_u_del = check_v3(p_u, SWE->p_u, SWE->NNodes);
      if (p_u_del < _TINY_VAL_){
	printf("...PASS\n");
      }
      else{
	printf("...FAIL\n");
	exit(-1);
      }

      printf("Computing Projection Operator (v) ");
      fType p_v_del = check_v3(p_v, SWE->p_v, SWE->NNodes);
      if (p_v_del < _TINY_VAL_){
	printf("...PASS\n");
      }
      else{
	printf("...FAIL\n");
	exit(-1);
      }

      printf("Computing Projection Operator (w) ");
      fType p_w_del = check_v3(p_w,SWE->p_w,SWE->NNodes);
      if (p_w_del < _TINY_VAL_){
	printf("...PASS\n");
      }
      else{
	printf("...FAIL\n");
	exit(-1);
      }
    }

  } 
  else{
    for(int i=0; i<3*SWE->NNodes; i++){
      SWE->p_u[i] = (fType) p_u[i];
      SWE->p_v[i] = (fType) p_v[i];
      SWE->p_w[i] = (fType) p_w[i];
    }
  }

  free(x);
  free(y);
  free(z);

  free(p_u);
  free(p_v);
  free(p_w);

  tmp = (double *) malloc(sizeof(double) * SWE->NNodes * 3);
  SWE->gradghm  = (fType*) malloc(sizeof(fType) * SWE->NNodes * 3);
  fread((void*) tmp, sizeof(double), SWE->NNodes * 3, fp_in);
  for(int i = 0; i < 3*SWE->NNodes; i++){
    SWE->gradghm[i] = (fType) tmp[i];
    }
  free(tmp);

  // allocate and read Derivative Matricess

  tmp    =  (double*) malloc(sizeof(double) * SWE->NNodes * SWE->NNbrs);
  SWE->Dx = (fType*) malloc(sizeof(fType) * SWE->NNodes * SWE->NNbrs_Padded);
  SWE->Dy = (fType*) malloc(sizeof(fType) * SWE->NNodes * SWE->NNbrs_Padded);
  SWE->Dz = (fType*) malloc(sizeof(fType) * SWE->NNodes * SWE->NNbrs_Padded);
  SWE->L =  (fType*) malloc(sizeof(fType) * SWE->NNodes * SWE->NNbrs_Padded);

  fread((void*) tmp, sizeof(double), SWE->NNodes * SWE->NNbrs, fp_in);
  for(int i=0; i<SWE->NNodes; i++){
    for(int j=0; j<SWE->NNbrs; j++){
      SWE->Dx[i*SWE->NNbrs_Padded + j] = (fType) tmp[i*SWE->NNbrs+j]/a;
    }
    for(int j=SWE->NNbrs; j<SWE->NNbrs_Padded; j++){
      SWE->Dx[i*SWE->NNbrs_Padded + j] = 0.0;
    }
  }   

  fread((void*) tmp, sizeof(double), SWE->NNodes * SWE->NNbrs, fp_in);
  for(int i=0; i<SWE->NNodes; i++){
    for(int j=0; j<SWE->NNbrs; j++){
      SWE->Dy[i*SWE->NNbrs_Padded + j] = (fType) tmp[i*SWE->NNbrs+j]/a;
    }
    for(int j=SWE->NNbrs; j<SWE->NNbrs_Padded; j++){
      SWE->Dy[i*SWE->NNbrs_Padded + j] = 0.0;
    }
  }   

  fread((void*) tmp, sizeof(double), SWE->NNodes * SWE->NNbrs, fp_in);
  for(int i=0; i<SWE->NNodes; i++){
    for(int j=0; j<SWE->NNbrs; j++){
      SWE->Dz[i*SWE->NNbrs_Padded + j] = (fType) tmp[i*SWE->NNbrs+j]/a;
    }
    for(int j=SWE->NNbrs; j<SWE->NNbrs_Padded; j++){
      SWE->Dz[i*SWE->NNbrs_Padded + j] = 0.0;
    }
  }   

  fread((void*) tmp, sizeof(double), SWE->NNodes * SWE->NNbrs, fp_in);
  for(int i=0; i<SWE->NNodes; i++){
    for(int j=0; j<SWE->NNbrs; j++){
      SWE->L[i*SWE->NNbrs_Padded + j] = (fType) tmp[i*SWE->NNbrs+j]*gamma;
    }
    for(int j=SWE->NNbrs; j<SWE->NNbrs_Padded; j++){
      SWE->L[i*SWE->NNbrs_Padded + j] = 0.0;
    }
  }   

  free(tmp);

   // read DM stencil ids (dimensions: NNodes x NNbrs)

  int *itmp = (int*) malloc(sizeof(int) * SWE->NNodes * SWE->NNbrs);
  SWE->idx = (int*) malloc(sizeof(int) * SWE->NNodes * SWE->NNbrs_Padded);
  fread((void*) itmp, sizeof(int), SWE->NNodes * SWE->NNbrs, fp_in);

  // 1-based indices in file must be shifted to 0-based indices

  for(int i=0; i<SWE->NNodes; i++){
    for(int j=0; j<SWE->NNbrs; j++){
      SWE->idx[i*SWE->NNbrs_Padded + j] = itmp[i*SWE->NNbrs + j] - 1;
    }
    for(int j=SWE->NNbrs; j<SWE->NNbrs_Padded; j++){
      SWE->idx[i*SWE->NNbrs_Padded + j] = 0;
    }
  }
  free(itmp);
  fclose(fp_in);

  // Read in initial conditions and correct values after 100 timesteps

  if (! strcmp(testcase,"tc5")){

    // read initial conditions from inputFile

    fType *Sfile = malloc( SWE->NState*SWE->NNodes * sizeof(fType) );
    SWE_read_initial_conditions(inputFile, Sfile);

    if (COMPUTE_TC5 == 1){

      // Compute initial conditions ab initio

      fType *Scomp = malloc( SWE->NState*SWE->NNodes * sizeof(fType) );

      SWE_init_testcase5(SWE, Scomp);

      fType S_del = check_state(Sfile, Scomp, SWE->NState, SWE->NNodes);
      printf("S_del = %.15f\n",S_del);

      if (S_del < _TINY_VAL_){

	// Copy Scomp into S

	for(int i=0; i<SWE->NState*SWE->NNodes; i++){
	  S[i] = Scomp[i];
	}
	free(Scomp);
	printf("...PASS\n");
      }
      else{
	printf("...FAIL\n");
	for(int i=0; i<SWE->NState*SWE->NNodes; i++){
	  S[i] = Sfile[i];
	}
	free(Sfile);
      }

    }
    else{

      // Copy Sfile into S

      for(int i=0; i<SWE->NState*SWE->NNodes; i++){
	S[i] = Sfile[i];
      }

      free(Sfile);

    }


  }
  else{
    printf("SWE_init ERROR: unsupported testcase %s, halting.\n",testcase);
    exit(-1);
  }
  printf("\nInitialization complete.\n");

}

