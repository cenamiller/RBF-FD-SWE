#include <stdlib.h>
#include <floating_types.h>
#include <profiling.h>
//#include <advance.h>

#ifdef CACHEQ
#include <cq.h>
#define getTime()         (cq_nstime() / 1e9)
#else
#define CQ_POOL(x)
#define CLANG_NOINLINE
#endif

extern timing_struct local_timer;

void halo_update(fType* S);

void RK3_advance(const fType dt, fType *t, const int NPts, fType *S, void f ( const fType t, const fType* S, fType* K)){

  const fType sixth = 1.0/6.0;

  fType *K1;
  fType *K2;
  fType *K3;

  fType *K = (fType*) malloc(NPts * sizeof(fType));
  fType *dS = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp1 = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp2 = (fType*) malloc(NPts * sizeof(fType));

  double t_start = getTime();
  K1 = K;
  f((*t), S, K1);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp1[i] =  S[i] + 0.5*dt*K1[i]; //2 flop/site
    Stmp2[i] =  S[i] - dt*K1[i]; //2 flop/site
    dS[i] = K1[i];
  }
  local_timer.t_update +=  (getTime() - t_start);
   
  t_start = getTime();
  K2 = K;
  f((*t)+0.5*dt, Stmp1, K2);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp2[i] += 2.0*dt*K2[i]; //2 flop/site
    dS[i] += 4.0*K1[i];  //2 flop/site
  }
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K3 = K;  
  f((*t)+dt, Stmp2, K3);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    S[i] += (sixth*dt)*(dS[i]+K3[i]); //2 flop/site
  }
  local_timer.t_update +=  (getTime() - t_start);

  (*t) += dt;

  free(K);
  free(dS);
  free(Stmp1);
  free(Stmp2);

}

#ifdef CACHEQ

//void RK4_advance(fType dt, CQ_timing_struct *t, int NPts, fType *S, void fswe( const CQ_timing_struct t, const fType* S, fType* K)){
void RK4_advance(fType dt, fType *t, int NPts, fType *S, void fswe( const fType t, const fType* S, fType* K)){

  const fType sixth = 1.0/6.0;

  fType *K1;
  fType *K2;
  fType *K3;
  fType *K4;

  fType *dS = (fType*) cq_malloc(3,NPts * sizeof(fType));
  fType *Stmp = (fType*) cq_malloc(3,NPts * sizeof(fType));
  fType *K = (fType*) cq_malloc(3,NPts * sizeof(fType));

  double t_start = getTime();
  K1 = K;
  fswe((*t),S,K1);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] =  S[i] + 0.5*dt*K1[i]; // 2 flop/site
    dS[i] = K1[i];
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K2 = K;
  fswe((*t)+0.5*dt, Stmp, K2);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + 0.5*dt*K2[i];  // 2 flop/site
    dS[i] += 2.0*K2[i]; // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K3 = K;
  fswe((*t)+0.5*dt, Stmp, K3);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + dt*K3[i];  // 2 flop/site
    dS[i] += 2.0*K3[i];         // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K4 = K;
  fswe((*t)+dt, Stmp, K4);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++)
    S[i] += (sixth*dt)*(dS[i]+K4[i]); //2 flop/site
  halo_update(S);
  local_timer.t_update +=  (getTime() - t_start);

  cq_free(K);
  cq_free(dS);
  cq_free(Stmp);

  (*t) += dt;
}

#else

void RK4_advance(const fType dt, fType *t, int NPts, fType *S, void f ( const fType t, const fType* S, fType* K)){

  const fType sixth = 1.0/6.0;

  fType *K1;
  fType *K2;
  fType *K3;
  fType *K4;

  fType *dS = (fType*) malloc(NPts * sizeof(fType));
  fType *Stmp = (fType*) malloc(NPts * sizeof(fType));
  fType *K = (fType*) malloc(NPts * sizeof(fType));

  double t_start = getTime();
  K1 = K;
  f((*t),S,K1);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] =  S[i] + 0.5*dt*K1[i]; // 2 flop/site
    dS[i] = K1[i];
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);
   
  t_start = getTime();
  K2 = K;
  f((*t)+0.5*dt, Stmp, K2);  
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + 0.5*dt*K2[i];  // 2 flop/site
    dS[i] += 2.0*K2[i]; // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K3 = K;
  f((*t)+0.5*dt, Stmp, K3);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++){
    Stmp[i] = S[i] + dt*K3[i];  // 2 flop/site
    dS[i] += 2.0*K3[i];         // 2 flop/site
  }
  halo_update(Stmp);
  local_timer.t_update +=  (getTime() - t_start);

  t_start = getTime();
  K4 = K;
  f((*t)+dt, Stmp, K4);
  local_timer.t_rhs +=  (getTime() - t_start);

  t_start = getTime();
  for(int i=0; i<NPts; i++)
    S[i] += (sixth*dt)*(dS[i]+K4[i]); //2 flop/site
  halo_update(S);
  local_timer.t_update +=  (getTime() - t_start);

  free(K);
  free(dS);
  free(Stmp);

  (*t) += dt;
}
#endif
