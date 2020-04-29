#!/bin/csh
source env.csh  
setenv SWE_NUM_COMPUTE_UNITS 1
setenv SWE_NODES 2562
setenv SWE_NSTEPS 100
setenv SWE_HIST_INT 1000
#setenv SWE_NATTEMPTS 10
setenv SWE_NATTEMPTS 1
setenv SWE_USE_RCM 1
setenv SWE_USE_NETCDF 0
setenv SWE_USE_OCL 0

setenv SWE_INPUT_FILE icos"$SWE_NODES"_tc5_input.bin

echo "Running gcc version of swe"
time ./pdex.gcc swe | tee out_gcc_2562.txt

echo "Running ncar31"
/opt/cacheq/qcc/0.5/bin/cqrun ncar31.cq swe | tee out_cq_2562.txt

echo "Running ncar31_sim"
/opt/cacheq/qcc/0.5/bin/cqrun ncar31_sim.cq swe -tp lx250 | tee out_cq_sim_2562.txt

