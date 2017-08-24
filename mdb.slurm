#!/bin/bash
#SBATCH -J mdb
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o mdb.o%j
#SBATCH -e mdb.e%j
#SBATCH -t 1:00:00
#SBATCH -A mega2014

module load launcher





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=mdb

$LAUNCHER_DIR/paramrun
