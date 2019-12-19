#!/bin/csh
   setenv SWE_NODES 2562
   if ("$1" != "") setenv SWE_NODES $1
   if ("$1" == 10) setenv SWE_NODES 10242
   if ("$1" == 40) setenv SWE_NODES 40962
   echo "Number of nodes is $SWE_NODES"
     
   setenv SWE_NUM_COMPUTE_UNITS 1
   setenv SWE_NSTEPS 100
   setenv SWE_HIST_INT 1000
   #setenv SWE_NATTEMPTS 10
   setenv SWE_NATTEMPTS 1
   setenv SWE_USE_RCM 1
   setenv SWE_USE_NETCDF 0
   setenv SWE_USE_OCL 0

   setenv SWE_INPUT_FILE $CACHEQ_SOURCE_DIR/customer/ncar.float32/input/icos"$SWE_NODES"_tc5_input.bin

   echo "===== Running gcc version of swe"
   time swe_gcc >& out.gcc

   echo "===== Running cq version of swe"
   setenv CQ_CACHE_MISS_DELAY 0
   setenv CQ_CACHE_ACCESSESPERSTRIPE 4
   setenv CQ_CACHE_NSTRIPES 32
   time $CACHEQ_SOURCE_DIR/cqrun/cqrun.py swe_cq |& tee out.cq
