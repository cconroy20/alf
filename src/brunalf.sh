#!/bin/bash

file='lris2.dat'
tag='test'

#get a list of the input data
arr=(`cat $file`)
len=${#arr[*]}

for i in `seq 1 $len`; do

cat >> tmp.${i}.slurm <<EOF
#!/bin/bash
#SBATCH -J alf_${arr[$i-1]}${tag}
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 48:00:00 # Runtime
#SBATCH -p conroy,general # Partition to submit to
#SBATCH --mem-per-cpu=1000 # Memory per cpu in MB (see also --mem)
#SBATCH -o /dev/null # Standard out goes to this file
#SBATCH -e /dev/null # Standard err goes to this file

cd /n/conroyfs1/cconroy/alf/src/
alf_varimf.exe ${arr[$i-1]} ${tag} > out.${arr[$i-1]}${tag}

EOF

sbatch tmp.${i}.slurm
rm -f tmp.${i}.slurm

done
