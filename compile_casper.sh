#!/bin/bash -l

#SBATCH -C v100
#SBATCH -A NCIS0002
#SBATCH -n 1
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=1g
#SBATCH -t 00:10:00

module load pgi/18.4

make EXEC=swe_default.exe CONFIG_DIR=pascal CFDL=0 SFDL=0 SPLIT_DEV=0
make EXEC=swe_sfdl.exe CONFIG_DIR=pascal CFDL=0 SFDL=1 SPLIT_DEV=0
make EXEC=swe_cfdl.exe CONFIG_DIR=pascal CFDL=1 SFDL=0 SPLIT_DEV=0
make EXEC=swe_cfdl_sfdl.exe CONFIG_DIR=pascal CFDL=1 SFDL=1 SPLIT_DEV=0
