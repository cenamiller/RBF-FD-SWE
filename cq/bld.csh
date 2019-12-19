#!/bin/csh
   set stripeopt = (-DSTRIPE_OPT)
   set ds = `echo "x$*" | fgrep 'DSTRIPE_OPT'`
   if ("$ds" != "") set stripeopt = (-DSTRIPE_OPT)
   set flags = (${stripeopt} -DCACHEQ -DUSE_BINIO -DUSE_CFDL -DPADDED_NUM=4 \
              -DSIMD_LENGTH=4 -DRUN_32BIT -Iinclude -I$CACHEQ_SOURCE_DIR/libcq )
   set gccfiles = (main.c profiling.c runtime_params.c input.c rcm.c  \
                   reorder_nodes.c layout.c matrix_transformations.c  \
                   init_patches.c halos.c)
   set ofiles = ""
   foreach f ($gccfiles)
      echo gcc -c -g $flags $f
      gcc -c -g $flags $f
      set ofiles = ($ofiles -l$f:r.o)
   end

   echo ===== Building GCC version =====
   echo gcc -c -g -o main_gcc.o -DRK_GCC_COMPILE $flags main.c
   gcc -c -g -o main_gcc.o -DRK_GCC_COMPILE $flags main.c
   echo gcc -c -g -O $flags rk4_rbffd_swe.c
   gcc -c -g -O $flags rk4_rbffd_swe.c
   set gppopts = (-o swe_gcc rk4_rbffd_swe.o \
                  `echo $gccfiles | sed 's/\.c/.o/g'` \
                  ${CACHEQ_BUILD_DIR}/libcq/libcq.a -lpthread)
   set gppopts = (`echo $gppopts | sed 's/main.o/main_gcc.o/'`)
   echo g++ $gppopts
   g++ $gppopts

   echo ===== Building QCC version =====
   set qccopts = (rk4_rbffd_swe.c -g guide.json $* -o swe_cq -e RK_substep \
                  $flags $ofiles --superNodesOff --qdir)
   echo qcc $qccopts
   $CACHEQ_SOURCE_DIR/qcc_py/qcc.py $qccopts
   echo 


