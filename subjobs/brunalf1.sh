#!/bin/bash
#SBATCH -J alf1
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 48:00:00 # Runtime
#SBATCH -p conroy,general # Partition to submit to
#SBATCH --mem-per-cpu=1000 # Memory per cpu in MB (see also --mem)
#SBATCH -o /dev/null # Standard out goes to this file
#SBATCH -e /dev/null # Standard err goes to this file

cd /n/conroyfs1/cconroy/alf/subjobs/
../bin/alf_single.exe n1407_rad02 test1 >& out.single

