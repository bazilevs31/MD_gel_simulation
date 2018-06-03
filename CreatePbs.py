#!/usr/bin/env python

# Program: CreatePbs.py
# Purpose: create pbs file to run on supercomputer
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreatePbs.py for help,
# example: python CreatePbs.py --run -name my_name_of_pbs
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


import os
import argparse
import termcolor

def read_pbs():
    """
    read information about the pbs file you want to create
    this read pbs can read what you want to do
    -about creating polymers
    -running the simulation including the equilibration
    -analazying the results, using certain program

    TODO: for analyzing part choose what kind of analysis you do want to do

    """
    parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # group = parser.add_mutually_exclusive_group(help='What do you want to do create melt or run sim')
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('--create', action='store_true',help='Create and equilibrate the run')
    # group.add_argument('--run', action='store_true', help='Run the simulation')
    # group.add_argument('--analyze', action='store_true', help='Analyze the alignment using AnalyzeChain.py')
    # group.add_argument('--cont', action='store_true', help='continue the simulation using in.continue')
    # group.add_argument('--stretch', action='store_true', help='do stretch')

    parser.add_argument("-p",
                        "--procs",
                        dest="procs",
                        default=12,
                        type=int,
                        help="Number of processors")
    parser.add_argument("-n",
                        "--nodes",
                        dest="nodes",
                        default=5,
                        type=int,
                        help="Number of nodes")
    parser.add_argument("-m",
                        "--memmory",
                        dest="memmory",
                        default=600,
                        type=int,
                        help="Memmory to use")
    parser.add_argument("-wm",
                        "--wallminutes",
                        dest="wallminutes",
                        default=0,
                        type=int,
                        help="Walltime minutes")
    parser.add_argument("-wh",
                        "--wallhours",
                        dest="wallhours",
                        default=8,
                        type=int,
                        help="Number of nodes")
    parser.add_argument("-name",
                        "--pbsname",
                        dest="pbsname",
                        default='poly',
                        type=str,
                        help="Name of pbs file")
    parser.add_argument("-ln",
                        "--lammpsfile",
                        dest="lammpsfile",
                        default='222gel_10nm_l43',
                        type=str,
                        help="Name of lammps file")
    parser.add_argument("-wn",
                        "--whole_nodes",
                        dest="whole_nodes",
                        default=False,
                        action='store_true',
                        help="to use or not the whole node")

    args = parser.parse_args()

    print args
    print termcolor.colored('parameters have been red', 'green')

    return args

def create_pbs(args):
    """
    produces stuff

    Parameters
    ----------
    args : dict
        argparse ditionary obtained by read_pbs

    Raises
    ------
    ValueError
        Description
    """
    lammpsexec = 'lmp_openmpi'
    mpirun = 'mpirun'
    import_libraries = 'module load library/openmpi/1.8.4-intel'
    loadscript = """\
    module load lammps/1Feb14
    module load intel/14.0.2
    """

    parameters = {'procs': args.procs,
                  'nodes': args.nodes,
                  'memmory': args.memmory,
                  'wallhours': args.wallhours,
                  'wallminutes': int(args.wallminutes),
                  'wallhours': int(args.wallhours),
                  'lammpsexec': lammpsexec,
                  'lammpsfile': args.lammpsfile,
                  'mpirun': mpirun,
                  'pbsname': args.pbsname,
                  'import_libraries':import_libraries}

    assert args.wallminutes <= 60, " wallminutes should be less than 60"

    heading_script = """\
#!/bin/bash
#PBS -S /bin/bash
#PBS -l pmem={0[memmory]}mb
""".format(parameters)
    if args.whole_nodes:
        heading_script += """#PBS -l procs={0[procs]}""".format(parameters)
    else:
        heading_script += """#PBS -l procs={0[procs]}""".format(parameters)
    heading_script += """\
#PBS -l walltime={0[wallhours]}:{0[wallminutes]}:00
#PBS -m bea
#PBS -M vasiliy.triandafilidi@gmail.com
cd $PBS_O_WORKDIR
""".format(parameters)
    runscript = """\
        module load library/openmpi/1.8.4-intel
        echo "starting " >> info.process
        {0[mpirun]} {0[lammpsexec]} -in {0[lammpsfile]} -echo log > log.{0[lammpsfile]}.txt
        RESULT=$?  # result of the execution of the latest command
        if [ $RESULT -eq 0 ]; then
          echo "successful simulation"
          rm *.pbs.*
        else
          echo failed
          return 1
        fi
        echo "finished simulation" >> info.process
    """.format(parameters)

    filename = "pbs_%s.pbs" % args.pbsname

    with open(filename, 'w') as lmp:
        lmp.write(heading_script + loadscript + runscript + '\n')


def main():
    """
    Reads parameters and runs create_pbs
    """
    print "I am doing reading paramters"
    args = read_pbs()
    create_pbs(args)

if __name__ == '__main__':
    main()
