#include <stdio.h>
#include <profiling.h>
#include <SWE.h>

void SWE_profiler(timing_struct local_timer, SWE_struct SWE, int RK_Order, int Nsteps){
  double rk_flop_per_step = RK_Order*((2*SWE.NNbrs-1)*4*(3+1) + 17*3 + 31)*SWE.NNodes;
  double t_rhs_per_step = local_timer.t_rhs/Nsteps;
  double rhs_Gflops = 1e-9*(rk_flop_per_step/t_rhs_per_step);
  double update_flop_per_step = (4+2*RK_Order)*(RK_Order*SWE.NNodes);
  double t_update_per_step = local_timer.t_update/Nsteps;
  double update_Gflops = 1.e-9*(update_flop_per_step/t_update_per_step);
  printf("\n==============================================================\n\n");
  printf("Benchmark results:\n\n");
  printf("time in rhs = %f Gflops = %f\n",t_rhs_per_step,rhs_Gflops);
  printf("time in update = %f Gflops = %f\n",t_update_per_step,update_Gflops);
  }
