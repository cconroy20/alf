#!/bin/bash
#SBATCH -J alf
#SBATCH -n 32 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 24:00:00 # Runtime
#SBATCH -p conroy # Partition to submit to
#SBATCH --mem-per-cpu=1000 # Memory per cpu in MB (see also --mem)
#SBATCH -o /dev/null # Standard out goes to this file
#SBATCH -e /dev/null # Standard err goes to this file

cd /n/conroyfs1/cconroy/alf/subjobs/
mpirun -n 32 ../bin/alf_simple_R5.exe megadeep Dec10_simple_R5 >& out.mega.simple_R5

