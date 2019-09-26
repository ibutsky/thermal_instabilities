#!/bin/bash
# Request three nodes that will not be allocated to others for the duration of the run.
#SBATCH --job-name=enzo_test
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=1:00:00


# Start from an "empty" module collection.
module purge
# Load in what we need to execute mpirun.
module load slurm gcc lib/hdf5 openmpi


cd /mnt/ceph/users/ibutsky/test/tcool_tff_1

c=0
while true
do
    if [ ! -e estd.out.$c ]; then
        break
    fi
    c=$[$c+1]
done
mv estd.out estd.out.$c


# Find out what nodes we were assigned.
#srun hostname

mpirun -n 16 ~/enzo-dev-cr/src/enzo/enzo.exe -d ThermalInstability.enzo >&estd.out
#mpirun -n 16 ~/enzo-dev-cr/src/enzo/enzo.exe -d -r DD0559/DD0559 >&estd.out  
