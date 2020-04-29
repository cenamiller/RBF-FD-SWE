// to build: /usr/bin/gcc main.c finish.c output.c advance.c rhs.c halo_update.c SWE_testcase5.c vmath.c Input.c profiling.c -O3 -I. -o pdex

#ifdef CACHEQ
#include "cq.h"   // use cq_malloc()
#define getTime() cq_nstime()
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <floating_types.h>
#include <SWE.h>
#include <NS.h>
#include <rcm.h>
#include <profiling.h>
#include <vmath.h>

// Only need _VERBOSE_ to debug initialization
#undef _VERBOSE_
#define NCheckAnswers 100

void RK3_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K));
#ifdef CACHEQ

//void RK4_advance(fType dt, CQ_timing_struct *t, int NPts, fType *S, void fswe( const CQ_timing_struct t, const fType* S, fType* K));
void RK4_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K));
#else
void RK4_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K));
#endif
void output(const fType t, const int NState, const fType *S);
void history(const char* pdecase, fType *S, fType *S_obs, int NPts);

fType* reorder_2D_fp_arr(fType* var, int* mapping, int dim1, int dim1_stride, int dim2);

void finish(fType *S, fType *K);

//supported RHSs

void fsimple( const fType t, const fType *S, fType* K );
void   faero( const fType t, const fType *S, fType* K );
#ifdef CACHEQ

//void    fswe( const CQ_timing_struct t, const fType* S, fType* K );
void    fswe( const fType t, const fType* S, fType* K );
#else
void    fswe( const fType t, const fType* S, fType* K );
#endif
void fstraka( const fType t, const fType* S, fType* K );

// SWE operators...

void SWE_init_bin(char* inputFile, char* testcase, SWE_struct *SWE, fType *S);
void SWE_get_init_filename(char *inputFile, int NNodes, char *testcase);
void SWE_profiler(timing_struct local_timer, SWE_struct SWE, int RK_Order, int Nsteps);
void SWE_correctness(SWE_struct SWE, char *inputFile, fType* S, fType* S_ref);
void SWE_final_conditions(char* inputFile, fType* S_ref);
  
// Global variables

SWE_struct* SWE;
NS_struct* NS;

timing_struct local_timer;

int main(int argc, char *argv[])
{

  int NState;
  int NNodes;
  int NPts;
  int Nsteps; 

  int RK_Order;
  fType dt;
  fType t;
  fType* S=NULL;
  fType* S_obs=NULL;
  fType* S_ref=NULL;
  fType* K=NULL;
  int *mapping=NULL;
  char inputFile[30];
  
  void (*fptr)( const fType t, const fType* S, fType* K) = NULL; //function pointer to assign RHS for specified problem
#ifdef CACHEQ
  CQ_timing_struct* cq_timing = (CQ_timing_struct*)malloc(sizeof(CQ_timing_struct));
  memcpy(cq_timing, &local_timer, sizeof(local_timer));
  cq_timing->t_update = cq_canonicalize_pointer(local_timer.t_update);
  cq_timing->t_rhs = cq_canonicalize_pointer(local_timer.t_rhs);
#endif
  if(argc==2)
    {
      if(! strcmp(argv[1], "simple")){
	NNodes = 1;    // number of nodes
	NState = 1;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100;  // number of timesteps
	RK_Order = 4;  // Runge Kutta Order
	dt = 0.1;      // timestep in seconds
	t = 0.0;
	S[0] = 1.0;
	fptr = &fsimple;
      }
      else if(! strcmp(argv[1], "aero")){
	NNodes = 1;    // number of nodes
	NState = 2;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100;  // number of timesteps
	RK_Order = 4;  // Runge Kutta Order
	dt = 0.1;      // timestep in seconds
	t = 0.0;
	S[0] = 2000.0;
	S[1] = 0.0;


	fptr = &faero;
      }
      else if(! strcmp(argv[1], "swe")){

	int RCM = 1;            // Reverse Cuthill McKee reorder 1 = Yes; 2 = No
	NNodes = 10242;         // Resolution in points
        t = 0;                  // Start time
	fptr = &fswe;       
	Nsteps = 200;           // number of timesteps
	RK_Order = 4;           // Runge Kutta Order

	// Build init file name

	char testcase[4];       
	strcpy(testcase,"tc5"); // Select testcase
	SWE_get_init_filename(inputFile, NNodes, testcase);
#ifdef _VERBOSE_
	printf("inputFile = %s\n", inputFile);
#endif
	NState = 4;
	NPts = NNodes*NState;

        SWE   = (SWE_struct*) malloc(sizeof(SWE_struct));
	S     = (fType*) malloc(NPts * sizeof(fType));
	S_obs = (fType*) calloc(NPts, sizeof(fType));
	S_ref = (fType*) calloc(NPts, sizeof(fType));

	K     = (fType*) malloc(NPts * sizeof(fType));
	SWE_init_bin(inputFile, testcase, SWE, S);
	dt = SWE->dt;           // timestep in seconds

	if (RCM == 1){
	  printf("\n==============================================================\
                  \n\nReordering RBF-FD stencil matrix using Reverse Cuthill McKee\
                  \n\nINITIAL stencil matrix:\n");
	  rcm_print_max_bandwidth(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);       
	  mapping = rcm_mapping(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);
	  rcm_check_mapping(mapping, SWE->NNodes);
	  SWE_reorder_nodes(SWE, mapping);
	  printf("\nREORDERED (RCM) stencil matrix:\n");
	  rcm_print_max_bandwidth(SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs, SWE->idx);       
	  S = reorder_2D_fp_arr(S, mapping, SWE->NNodes, NState, NState);
	}
	else{
	  for(int i=0; i<NNodes; i++){
	    mapping[i]=i;
	  }
	}

	// Initialize verification array for SWE case

	SWE_final_conditions(inputFile, S_ref);
	printf("\nReordering verification State Vector");
	S_ref = reorder_2D_fp_arr(S_ref, mapping, SWE->NNodes, NState, NState);
	printf(" ...done\n");
#ifdef _VERBOSE_
	output(t,NState,S);
#endif
      }
      else if(! strcmp(argv[1], "straka")){
        NS= (NS_struct*) malloc(sizeof(NS_struct));
	NNodes = 3000; // number of RBF nodes
	NState = 4;    // number of state variables in SW PDE
	NPts   = NState*NNodes;
	S = (fType*) malloc(NPts * sizeof(fType));
	Nsteps = 100; // number of timesteps
	RK_Order = 4; // Runge Kutta Order
	dt = 2;       // timestep in seconds
	t = 0.0;
	fptr = &fstraka;
      }
      else{
	printf("ERROR IN MAIN: Unknown PDE case %s\n",argv[1]);
	exit(-1);
      }
    }
  else{
    printf("ERROR IN MAIN: Invalid num of command line arguments\n");
    printf("Correct usage is: ./pdex {simple,aero,swe,straka}\n");    
    exit(-1);
  }

  init_timer(&local_timer);

#ifdef _VERBOSE_
  output(t,NState,S);
#endif

  // Integrate the equations...
  printf("\n==============================================================\n\n");
  printf("Advancing PDE equations %d timesteps with Runge Kutta Order = %d\n",Nsteps,RK_Order);

  if (RK_Order == 3){
    for(int it=0; it<Nsteps; it++){
      RK3_advance(dt,&t,NPts,S,fptr);
      if (it == NCheckAnswers-1){
	printf("\nSaving state variable After %d iterations ",it+1);
	history(argv[1],S,S_obs,NPts);
	printf("...done\n");
      }
#ifdef _VERBOSE_
      output(t,NState,S);
#endif
    }
  }
  else if (RK_Order == 4){
    for(int it=0; it<Nsteps; it++){
#ifdef CACHEQ

      //RK4_advance(dt,cq_timing,NPts,S,fswe); //S-Size of NPts
      RK4_advance(dt,&t,NPts,S,fswe); //S-Size of NPts
#else
      RK4_advance(dt,&t,NPts,S,fswe); //S-Size of NPts
#endif
      if (it==NCheckAnswers-1){
	printf("\nSaving state variable After %d iterations ",it+1);
	history(argv[1],S,S_obs,NPts);
	printf("...done\n");
      }
#ifdef _VERBOSE_
      output(t,NState,S);
#endif
    }
  }
  else{
    printf("ERROR IN MAIN: RK_Order = %i not supported\n",RK_Order);
    exit(-1);
  }

  if(! strcmp(argv[1], "swe")){
    SWE_profiler(local_timer,*SWE,RK_Order,Nsteps);
    SWE_correctness(*SWE, inputFile, S_obs, S_ref);
    }     

  finish(S,S_ref);
}
