#!/bin/bash
#SBATCH -n 6
#SBATCH --ntasks-per-node=6
#SBATCH -c 4
#SBATCH -t 01:00:00

module load cesga/2020
module load intel impi

export HOME=/home/ulc/es/pgg/aco-cellnopt
export NUM_PROC_MPI=6
export OMP_NUM_THREADS=4

srun $HOME/build/aco $HOME/benchmark/LiverDREAM/insilico_model_data.h5


