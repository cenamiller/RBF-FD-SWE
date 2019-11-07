#!/bin/bash -l

#SBATCH -C v100
#SBATCH -A NCIS0002
#SBATCH -n 4
#SBATCH --exclusive
#SBATCH --gres=gpu:v100:4
#SBATCH --mem=1g
#SBATCH -t 00:10:00

export KMP_AFFINITY=disabled  #granularity=core,balanced
export OMP_SCHEDULE=static
export OMP_NUM_THREADS=1

export SWE_NODES=40962
export SWE_NSTEPS=100
export SWE_HIST_INT=1000
export SWE_NATTEMPTS=1
export SWE_USE_RCM=1
export SWE_USE_NETCDF=0
export SWE_USE_OCL=0
export SWE_NUM_COMPUTE_UNITS=1
export SWE_OUTPUT_NAME="output"

export SWE_TOP_DIR=$PWD
export SWE_INPUT_FILE=$SWE_TOP_DIR/input_files/icos"$SWE_NODES"_tc5_input.bin
export SWE_RUN_DIR=$SWE_TOP_DIR/run/psg/

module load pgi/18.4
module load openmpi

for Layout in cfdl sfdl cfdl_sfdl default; do
export NumSubdevices=1
while [ $NumSubdevices -le 4 ]; do
echo >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt
echo "Run script name: run.sh" >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt
echo "NumSubdevices = $NumSubdevices" >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt

if [ $NumSubdevices -gt 30 ]; then
export SWE_NATTEMPTS=30
elif [ $NumSubdevices -gt 16 ]; then
export SWE_NATTEMPTS=20
else
export SWE_NATTEMPTS=10
fi

export num_tasks=$NumSubdevices
mpirun -np $NumSubdevices $SWE_RUN_DIR/gpu_launch.sh $SWE_RUN_DIR/swe_"$Layout".exe >> "$SWE_RUN_DIR"/output/"$SWE_OUTPUT_NAME"_"$Layout".txt

export NumSubdevices=$[$NumSubdevices+1]
done
done
echo "DONE!"
