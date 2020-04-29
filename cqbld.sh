#!/bin/csh
make clean
source env.csh   
   #set cfiles = (main.c finish.c output.c advance.c rhs.c halo_update.c SWE_init_testcase5.c vmath.c Input.c profiling.c) 
   set cfiles = (finish.c fstraka.c halo_update.c history.c main.c output.c profiling.c rcm.c reorder_nodes.c rhs.c  SWE_correctness.c SWE_final_conditions.c SWE_get_init_filename.c SWE_init_bin.c SWE_init_testcase5.c SWE_profiler.c SWE_read_initial_conditions.c SWE_reorder_nodes.c vmath.c)
   echo "====== Building GCC version ======"
   set ofiles = 
   foreach f (advance.c $cfiles)
      echo gcc -c -DUSE_CFDL -O $f -I include -I/opt/cacheq/qcc/0.5/
      gcc -c -DUSE_CFDL -O $f -I include -I/opt/cacheq/qcc/0.5/ -Wno-unused-result
      set ofiles = ($ofiles $f:r.o)
   end
   g++ -o pdex.gcc $ofiles /opt/cacheq/qcc/0.5/lib/libcq.a -lpthread -lm
   
   echo 
   echo "====== Building CacheQ version ======"
   set ofiles = 
   foreach f ($cfiles)
      echo gcc -DUSE_CFDL -c $f -I include -DCACHEQ -I/opt/cacheq/qcc/0.5/
      gcc -DUSE_CFDL -c $f -I include -DCACHEQ -I/opt/cacheq/qcc/0.5/

      set ofiles = ($ofiles -l $f:r.o)
   end

   set cmdline = (/opt/cacheq/qcc/0.5/bin/qcc advance.c -I include -D USE_CFDL -o ncar31 -e RK4_advance $ofiles $*)
   echo $cmdline
   $cmdline
   echo "===== Building Simulation version ======"

   set cmdline = (/opt/cacheq/qcc/0.5/bin/qcc advance.c -t -I include -D USE_CFDL -o ncar31_sim -e RK4_advance $ofiles $*)
   echo $cmdline
   $cmdline

   
