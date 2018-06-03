
   lmp_mpi -in 1_equil_lj.lammps -var infile filename -var Temp 1.1 -var kspring 30  -var name t11_ljeq1 -var nsteps 10000 -echo both

   lmp_mpi -in 2_equil_coul.lammps -var infile filename_t11_ljeq1 -var Temp 1.1 -var kspring 30  -var name coul_eq  -var nsteps 10000 -echo both


   lmp_mpi -in new_coul_run.lammps -var infile filename_t11_ljeq1_coul_eq -var Temp 1.1 -var kspring 30  -var name coulsim -var nsteps 10000  -echo both

  lmp_mpi -in pressure_profile.lammps -var infile filename_t11_ljeq1_coul_eq_coulsim -var Temp 1.1 -var kspring 30  -var name press -echo both
