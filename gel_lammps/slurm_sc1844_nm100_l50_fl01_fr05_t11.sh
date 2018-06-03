#!/bin/bash
#SBATCH --ntasks 121               # number of MPI processes
#SBATCH --cpus-per-task 1        # number of OpenMP threads per MPI process
#SBATCH --mem-per-cpu 600       # memory limit per CPU core (megabytes)
#SBATCH --time 6:50:00           # time limit (D-HH:MM:ss)
#SBATCH --mail-user=vasiliy.triandafilidi@gmail.com
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
module load intel/2016.4 openmpi/2.1.1
module load lammps-user-intel/20170331
#seems like 5 hours per 1M steps is enough
#so 15 hours per 3M steps should be enough

 

 srun  lmp_intel_cpu_openmpi -in 1_equil_lj.lammps -var infile sc1844_nm100_l50_fl01_fr05 -var Temp 1.1 -var kspring 30  -var name t11_ljeq1 -var nsteps 10000 -echo both

 srun  lmp_intel_cpu_openmpi -in 2_equil_coul.lammps -var infile sc1844_nm100_l50_fl01_fr05_t11_ljeq1 -var Temp 1.1 -var kspring 30  -var name coul_eq  -var nsteps 10000 -echo both


 srun  lmp_intel_cpu_openmpi -in new_coul_run.lammps -var infile sc1844_nm100_l50_fl01_fr05_t11_ljeq1_coul_eq -var Temp 1.1 -var kspring 30  -var name coulsim -var nsteps 10000  -echo both

 srun lmp_intel_cpu_openmpi -in pressure_profile.lammps -var infile sc1844_nm100_l50_fl01_fr05_t11_ljeq1_coul_eq_coulsim -var Temp 1.1 -var kspring 30  -var name press -echo both

#srun lmp_intel_cpu_openmpi -in pressure_profile.lammps -var infile sc1844_nm100_l50_fl01_fr05_t11_ljeq1_coul_eq_coulsim -var Temp 1.1 -var kspring 30  -var name pressv2 -echo both

#---- 
