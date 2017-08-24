#!/bin/bash
#SBATCH -J blast
#SBATCH -n 190
#SBATCH -N 16
#SBATCH -p normal
#SBATCH -o blast.o%j
#SBATCH -e blast.e%j
#SBATCH -t 3:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mstudiva@fau.edu

module load launcher





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=bl

$LAUNCHER_DIR/paramrun
