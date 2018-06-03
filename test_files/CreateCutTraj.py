#!/usr/bin/env python

# Program: CreateCutTraj.py
# Purpose: Cuts Trajectory for download (in future will wrap it as well)
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreateCutTraj.py -h for help, 
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 

# todo: add wrapping feature, make function look nice

from read_parameters import read_traj_vmd
import os
from termcolor import colored


args = read_traj_vmd()
psffile = os.path.splitext(args.psffile)[0]
PSF = psffile+'.psf'
DCD = os.path.splitext(args.traj)[0] + '.dcd'
DATA = os.path.splitext(args.datafile)[0] + '.data'
vmdfile = "create_cuttraj.vmd"
bashfile = "create_cuttraj.sh"

resfolder = psffile+'_results'
# write VMD loader script
parameters = {'vmdfile': vmdfile,
              'topology': PSF,
              'datafile': DATA,
              'trajectory': os.path.basename(DCD),
              'trajskip':args.trajskip,
              'resfolder':resfolder,
              'bashfile' : bashfile}

script = """\
mol new "{0[topology]}"
animate read dcd "{0[trajectory]}" waitfor all
animate write dcd skip{0[trajectory]} skip 200
echo "writing dcd to skip{0[trajectory]} "
exit
""".format(parameters)

bashscript = """\
mkdir {0[resfolder]}
vmd -dispdev text -e {0[vmdfile]}
cp skip{0[trajectory]} ./{0[resfolder]}
cp ./figures/*.pdf ./{0[resfolder]}
cp ./figures/*.npz ./{0[resfolder]}
# cp ./figures/log* ./{0[resfolder]}
cp *.psf ./{0[resfolder]}
""".format(parameters)

with open(vmdfile, 'w') as tcl:
    tcl.write(script+'\n')    
with open(bashfile, 'w') as tcl:
    tcl.write(bashscript+'\n')

os.system("bash {0[bashfile]}".format(parameters))
print colored("If there is an error with {0}: 'source {0}' to load everything manually, then repeat ".format(vmdfile),"blue")
print colored("If there is an error with {0}: 'bash {0}' to load everything manually, then repeat ".format(bashfile),"blue")
print colored("running the python script with explicict parameters that were generated".format(vmdfile),"blue")