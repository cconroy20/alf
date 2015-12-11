#!/bin/bash

file='inlist_n1600.dat'
tag='Aug05'
tag2=''

#get a list of the input data
arr=(`cat $file`)
len=${#arr[*]}
$len=1

for i in `seq 1 $len`; do

cat >> tmp.${i}.slurm <<EOF
#!/bin/bash
#SBATCH -J alf_${arr[$i-1]}${tag}
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 24:00:00 # Runtime
#SBATCH -p conroy # Partition to submit to
#SBATCH --mem-per-cpu=1000 # Memory per cpu in MB (see also --mem)
#SBATCH -o /dev/null # Standard out goes to this file
#SBATCH -e /dev/null # Standard err goes to this file

cd /n/conroyfs1/cconroy/alf/subjobs/
mpirun -n 16 ../bin/alf.exe ${arr[$i-1]}${tag2} ${tag} >& out.${arr[$i-1]}${tag}

EOF

sbatch tmp.${i}.slurm
rm -f tmp.${i}.slurm

done
