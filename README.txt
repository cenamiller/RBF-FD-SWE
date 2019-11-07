
------------------------------------------------------------------------------------------------------------------

Radial Basis Functions with Finite Differencing for Shallow Water Equations

Last update: Richelle Streater, September 7, 2018

------------------------------------------------------------------------------------------------------------------

Getting started

Once all folders are on device, navigate to top directory: Gen3/

On HPCL:
Compile with the following command: sbatch compile_hpcl.sh
Run with the following command: sbatch run/HPCL/runCL.sh from Gen3/ directory
Expected output: output_cfdl.txt, output_sfdl.txt, output_cfdl_sfdl.txt, output_default.txt in run/HPCL/output/

On personal device:
Compile with the following command: . ./compile.sh
Run with the following command: . ./run/openCL/run.sh
Expected output: output_cfdl.txt, output_sfdl.txt, output_cfdl_sfdl.txt, output_default.txt in run/openCL/output/

------------------------------------------------------------------------------------------------------------------

Requirements

To run with OpenCL:
--> OpenCL Version 2.0 or higher
--> Add -L and -I flags to OCL_LIBS/OCL_FLAGS to config.swe if necessary
--> Can set SPLIT_DEV=0 if device does not support subdividing
--> Set OPENCL=1 in config.swe before compiling
--> Set SWE_USE_OCL=1 in run script

To run with OpenMP:
--> Set OPENMP=1 in config.swe before compiling
--> Set OMP_NUM_THREADS in run script

To run with MPI:
--> Set MPI=1 in config.swe before compiling
--> OpenCL with multiple tasks is possible, but device splitting with MPI is not implemented
--> Change LD_LIBRARY_PATH and PATH in run/hpcl/runMP.sh or run/hpcl/runCL.sh if necessary

To use NetCDF:
--> Set NCIO=1 in config.swe before compiling
--> Set SWE_USE_NETCDF=1 in run script and set SWE_INPUT_FILE to a .nc file
--> Change NETCDF variable in include.mk if necessary

To run with Intel Compiler:
--> Set MPICC=mpiicc and CC=icc in config.swe
--> Load icc module before compiling if necessary
--> Change LD_LIBRARY_PATH and PATH in run/hpcl/runMP.sh or run/hpcl/runCL.sh if necessary

------------------------------------------------------------------------------------------------------------------

Replicating test results in read_output_file/results.xlsx on NCAR HPCL

For all:
--> To compile: "sbatch compile_hpcl.sh"
--> To run: "sbatch run/hpcl/runMP.sh" or "sbatch run/hpcl/runCL.sh" from Gen3/ directory
--> To get output files: Copy output files into read_output_file/output and run read_output_file/read_script.cpp

Tabs 10242 through 655362:
--> Set OPT_FLAGS = -O3 in arch/hpcl/config.swe
--> Compile with "CC=gcc" and "MPICC=mpicc" in arch/hpcl/config.swe
--> In run/HPCL/runCL.sh, set SWE_NODES to desired number (ex. 10242)
--> In run/HPCL/runCL.sh, Layout array should be "cfdl sfdl cfdl_sfdl default" and SWE_NODES=40962
--> Run with runCL.sh

OpenMP gcc tab:
--> Set OPT_FLAGS = -O3 in arch/hpcl/config.swe
--> Compile with "CC=gcc" and "MPICC=mpicc" in arch/hpcl/config.swe
--> In run/HPCL/runMP.sh, Layout array should be "cfdl sfdl cfdl_sfdl default" and SWE_NODES=40962
--> Run

OpenMP icc tab:
--> Set OPT_FLAGS = -O3 -xHost in arch/hpcl/config.swe
--> Compile with "CC=icc" and "MPICC=mpiicc" in arch/hpcl/config.swe
--> Repeat steps 1-6 in "OpenMP gcc tab"

OpenMP KMP_aff tab:
--> Compile with "CC=icc" and "MPICC=mpiicc" in arch/hpcl/config.swe
--> In run/HPCL/runMP.sh, set Layout array to "cfdl" and SWE_NODES=40962
--> Set KMP_AFFINITY=compact and run
--> Set KMP_AFFINITY=disabled and run
--> Set KMP_AFFINITY=scatter and run
--> Set KMP_AFFINITY=balanced and run

Aliasing tab:
--> Set "CC=gcc" and "MPICC=mpicc" in arch/hpcl/config.swe
--> Compile with OPT_FLAGS=-fstrict-aliasing in arch/hpcl/config.swe
--> In run/HPCL/runMP.sh, set Layout array to "cfdl" and SWE_NODES=40962
--> Run
--> Compile with OPT_FLAGS=-fno-strict-aliasing in arch/hpcl/config.swe and run with runMP.sh

------------------------------------------------------------------------------------------------------------------

Directory structure

Top Directory Folders:

Gen3/arch: Contains configuration parameters for hpcl and pascal testing and for general 
gnu or intel setup

Gen3/inputFiles: Contains all binary/netcdf input files for code

Gen3/read_output_file: Contains code to read eval_rhs values from output files

Gen3/run: Contains run scripts for HPCL and general OpenCL setup

Gen3/swe_code: Contains all c/cl code

------------------------------------------------------------------------------------------------------------------

swe_code folder structure:

Gen3/swe_code/io:
--> input.c:  Reads input files, either with binary or NetCDF format, and fills all differentiation matrices, 
state variable matrices, ordering, and constants
--> nc2bin.c: Converts to binary file from .nc format (not called by main function)

Gen3/swe_code/layout:
--> layout.c: Calls padding/reordering functions for differentiation matrices/state variable matrices
--> matrix_transformations.c: Functions to pad matrices (to allow for tiling/vectorizing) and rearrange based 
on CFDL/SFDL options

Gen3/swe_code/main:
--> main.c: Calls reading/reordering functions and calls patch initialization functions (for MPI). Declares
OpenCL objects and compiles kernels, opens device/platform, loads buffers, and sets kernel arguments for 
OpenCL. For n attempts and time steps, calls Runge-Kutta stepping function. Compares results to known array.
--> profiling.c: Use arrays of loop times to determine average time, min/max, and std dev for all operations.
Prints timing results.
--> rk4_rbffd_swe.c: Computes Runge-Kutta step with radial basis function finite differencing algorithm.
--> runtime_params.c: processes external variables set in run script.

Gen3/swe_code/mpi:
--> halos.c: Function for exchanging neighbor node information so that state variable matrix can be divided
among MPI threads
--> init_patches.c: Creates divided matrices and copies read and reordered arrays from thread 0 to other 
MPI threads.

Gen3/swe_code/ocl:
--> buffers.c: converts arrays into OpenCL buffer objects to be passed into the kernels and frees buffers
--> device_setup.c: Creates all OpenCL objects: kernels, devices, platforms, and command queues
--> RK_ocl.c: Version of rk4_rbffd_swe.c that used OpenCL and calls openCL kernels
--> kernel.cl: All Runge-Kutta step functions; vectorized along nodes and works for all layouts
--> kernel_CFDL.cl: All RK step functions; vectorized along u,v,w,h in state variable matrices. Only valid for
CFDL layout.

Gen3/swe_code/rcm:
--> rcm.c: Calls all reorder functions for Reverse CutHill-McKee ordering scheme
--> reorder_nodes: Defines mapping for Reverse CutHill-McKee ordering scheme

------------------------------------------------------------------------------------------------------------------