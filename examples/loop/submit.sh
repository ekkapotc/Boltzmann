#!/usr/bin/env zsh

#BSUB -J LOOP.N200.L5000.SERIAL_BLOCK_BLKSIZE_22ITERS
#BSUB -o LOOP.%J
#BSUB -W 00:30
#BSUB -M 1024
#BSUB -N
#BSUB -a openmpi
#BSUB -n 1
#BSUB -m mpi-s
#BSUB -R "span[ptile=1]"
#BSUB -x

time $MPIEXEC $FLAGS_MPI_BATCH ./main 200 5000 22
