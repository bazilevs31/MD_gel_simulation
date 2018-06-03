#!/bin/bash
#SBATCH --ntasks {0[procs]}               # number of MPI processes
#SBATCH --cpus-per-task 1        # number of OpenMP threads per MPI process
#SBATCH {0[mem_string]}
#SBATCH --time {0[wallhours]}:{0[wallminutes]}:00           # time limit (D-HH:MM:ss)
#SBATCH --mail-user=vasiliy.triandafilidi@gmail.com
#SBATCH --mail-type=ALL
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-1}}"
