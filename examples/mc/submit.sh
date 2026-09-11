#!/usr/bin/env zsh

#BSUB -J  MC.N1000.M2000.ITER4
#BSUB -o MC.N1000.M2000.ITER4%J
#BSUB -W  06:00
#BSUB -M  5120
#BSUB -N
#BSUB -n  1
#BSUB -m mpi-bull
#BSUB -a openmpi

time $MPIEXEC $FLAGS_MPI_BATCH ./main scenario_6.in 1000 2000 848454720
