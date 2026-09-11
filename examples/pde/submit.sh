#!/usr/bin/env zsh

#BSUB -P rwth0166
#BSUB -J PDE.M1000.N1000
#BSUB -o PDE.M1000.N1000.%J
#BSUB -W 01:00
#BSUB -M 1024
#BSUB -N
#BSUB -n 502
#BSUB -a openmpi
#BSUB -m mpi-bull

time $MPIEXEC $FLAGS_MPI_BATCH ./main scenario_6.in 1000 1000 
