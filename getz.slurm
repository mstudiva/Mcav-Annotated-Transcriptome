#!/bin/bash
#SBATCH -J getz
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o getz.o%j
#SBATCH -e getz.e%j
#SBATCH -t 2:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mstudiva@fau.edu

module load launcher





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=getz

$LAUNCHER_DIR/paramrun
