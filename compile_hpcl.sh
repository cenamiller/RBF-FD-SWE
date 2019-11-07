#!/bin/bash

export LM_LICENSE_FILE="28518@128.117.177.41"
export INTEL_LICENSE_FILE="28518@128.117.177.41"

source /usr/local/intel/2018u3/bin/compilervars.sh intel64

make EXEC=swe_default.exe CONFIG_DIR=hpcl CFDL=0 SFDL=0 RHS_SIMD_METHOD=2 GCC=0
make EXEC=swe_sfdl.exe CONFIG_DIR=hpcl CFDL=0 SFDL=1 RHS_SIMD_METHOD=2 GCC=0
make EXEC=swe_cfdl.exe CONFIG_DIR=hpcl CFDL=1 SFDL=0 RHS_SIMD_METHOD=2 GCC=0
make EXEC=swe_cfdl_sfdl.exe CONFIG_DIR=hpcl CFDL=1 SFDL=1 RHS_SIMD_METHOD=2 GCC=0

mpicc --showme:command
