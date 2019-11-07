#!/bin/bash

#SBATCH -n 1
#SBATCH -c 54
#SBATCH --wait
#SBATCH --mem-per-cpu=2000
#SBATCH --exclusive
#SBATCH -p skylake
#SBATCH -t 00:20:00

export SWE_NODES=40962
export SWE_NSTEPS=100
export SWE_HIST_INT=1000
export SWE_NATTEMPTS=20
export SWE_USE_RCM=1
export SWE_USE_NETCDF=0
export SWE_USE_OCL=1
export SWE_NUM_COMPUTE_UNITS=54
export SWE_OUTPUT_NAME="output"

export SWE_TOP_DIR=$PWD
export SWE_INPUT_FILE=$SWE_TOP_DIR/input_files/icos"$SWE_NODES"_tc5_input.bin
export SWE_RUN_DIR=$SWE_TOP_DIR/run/HPCL/

export LD_LIBRARY_PATH=/usr/local/opencl/ocl-icd/2.2.12/lib:/usr/local/opencl/icd/intel/cpu/6.4.0.37/opt/intel/opencl-1.2-6.4.0.37/lib64:/opt/slurm/latest/lib64:/usr/local/openmpi/3.1.0/lib:$LD_LIBRARY_PATH:$LD_LIBRARY_PATH
cd $SWE_TOP_DIR
source /usr/local/intel/2018u3/bin/compilervars.sh intel64
export PATH=/opt/slurm/latest/bin:/usr/local/openmpi/3.1.0/bin/:$PATH

for Layout in cfdl sfdl cfdl_sfdl default; do
export NumSubdevices=2
while [ $NumSubdevices -le 4 ]; do

echo >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt
echo "Run script name: runCL.sh" >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt

if [ $NumSubdevices -gt 30 ]; then
export SWE_NATTEMPTS=30
elif [ $NumSubdevices -gt 16 ]; then
export SWE_NATTEMPTS=20
else
export SWE_NATTEMPTS=10
fi

export SWE_NUM_COMPUTE_UNITS=$NumSubdevices
srun -o --mpi=pmi2 -c $NumSubdevices $SWE_RUN_DIR/swe_"$Layout".exe >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt

export NumSubdevices=$[$NumSubdevices+2]
done
done
echo "DONE!"
