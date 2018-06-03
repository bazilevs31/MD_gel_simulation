#!/usr/bin/env python

# Program: Create_VMD_dxanalyze_v2.py
# Purpose: creates a vmd script to analyze the system
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreCreate_VMD_analyze_dx.py -h for help,

# Copyright (c) 2018 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

from gridData import Grid
import os.path
import termcolor
import os
import argparse



# def create_slurm_string(procs, wh, dataname, tempname, temp):
def create_slurm_string(**params):

    runscript = """ """
    runscript += """\
    vmd -dispdev text -e {0[vmd_out_name]}.vmd
#---- """.format(params)

    heading_script = """\
#!/bin/bash
#SBATCH --ntasks {0[procs]}               # number of MPI processes
#SBATCH --cpus-per-task 1        # number of OpenMP threads per MPI process
#SBATCH --mem=82G
#SBATCH --time {0[wallhours]}:{0[wallminutes]}:00           # time limit (D-HH:MM:ss)
#SBATCH --mail-user=vasiliy.triandafilidi@gmail.com
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-1}}"
module load vmd/1.9.3

""".format(params)

    return heading_script  + runscript

def write_file(**params):
    """writes the vmd file"""


    if not os.path.exists(params['out_folder']):
        os.makedirs(params['out_folder'])

    headstring = """
    package require pbctools
    package require topotools"""


    body_string = """
    topo readlammpsdata {0[datafile]} full waitfor all

    animate read dcd {0[dcdfile]} waitfor all
    package require pmepot
    pbc wrap
    # # pmepot -mol 0  -ewaldfactor 0.25
    pmepot -mol {0[molnumber]} -ewaldfactor {0[ewaldfactor]}  -frames all  -grid {{ {0[nx]} {0[ny]} {0[nz]} }}    -dxfile {0[out_folder]}/{0[dxname]}.dx -loadmol none
    mol delete {0[molnumber]}

    """.format(params)


    with open("{0[vmd_out_name]}.vmd".format(params), 'a') as f:
        f.write(headstring  + body_string +  '\n' +'\n')
    return None




def run_analysis(args):
    """runs the actual analysis"""

    params_dict = vars(args)


    with open("{0[vmd_out_name]}.vmd".format(params_dict), 'w') as f:
        f.write('\n')

    temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
    molnumber = 0
    # for filename, _ in file_marker_dict.items():
    for filename in params_dict['files']:
        for tempname in temp_name_array:
            params_dict['filename'] = filename
            params_dict['tempname'] = tempname
            dcdfile = '{0[dcd_folder]}/{0[filename]}_t{0[tempname]}_{0[dcd_suffix]}.dcd'.format(params_dict)
            datafile = '{0[data_folder]}/{0[filename]}_t{0[tempname]}_{0[data_suffix]}.data'.format(params_dict)
            dxname = '{0[filename]}_t{0[tempname]}'.format(params_dict)

            if os.path.isfile(dcdfile) and os.path.isfile(datafile):
                print('file exists')
                params_dict['datafile'] = datafile
                params_dict['dcdfile'] = dcdfile
                params_dict['dxname'] = dxname
                params_dict['molnumber'] = molnumber
                print(termcolor.colored('adding VMD entry with the parameters', 'blue'))
                print(params_dict)

                write_file(**params_dict)
                molnumber += 1

        # finishing the writing with the exit
    with open("{0[vmd_out_name]}.vmd".format(params_dict), 'a') as f:
        f.write(' \n exit \n')

    if params_dict['create_job_flag']:
        with open("slurm_{0[vmd_out_name]}.sh".format(params_dict), 'w') as f:
            f.write(create_slurm_string(**params_dict))
            f.write(' \n exit \n')

    return None


if __name__ == '__main__':
    def is_valid_file(parser, arg):
        """
        Check if arg is a valid file
        """
        arg = os.path.abspath(arg)
        if not os.path.exists(arg):
            parser.error("The file %s doesn't exist " % arg)
        else:
            return arg


    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # datadir_to_walk = './data/'
    # dcd_filter = 'sc*press.dcd'
    # dcd_suffix = '_ljeq1_coul_eq_coulsim_press'
    # _vmdfile = 'glob_vmd_all'

    parser.add_argument('-a', action='append', dest='files',
                        type=str,
                        default=[],
                        help='Add values to analyze values to a list type (default: %(default)s)',
                        )

    parser.add_argument('-od', "--out_folder", dest="out_folder",
                            default='./dx',
                            help="where to output the dx files trajectory files directory)")

    parser.add_argument('-d', "--dcd_folder", dest="dcd_folder",
                            default='./data',
                            help="Input trajectory files directory)")

    parser.add_argument('-dd', "--data_folder", dest="data_folder",
                            default='./data',
                            help="Input trajectory files directory)")

    parser.add_argument("-ds",
                        "--dcd_suffix",
                        dest="dcd_suffix",
                        default='ljeq1_coul_eq_coulsim_press',
                        type=str,
                        help="suffix to cut from dcd file (default: %(default)s)")

    parser.add_argument("-dats",
                        "--data_suffix",
                        dest="data_suffix",
                        default='ljeq1_coul_eq_coulsim',
                        type=str,
                        help="suffix to cut from data file (default: %(default)s)")

    parser.add_argument("-o",
                        "--vmd_out_name",
                        dest="vmd_out_name",
                        default='all_vmd_files',
                        type=str,
                        help="the name of output vmd (default: %(default)s)")

    parser.add_argument("-ef", "--ewaldfactor", dest="ewaldfactor",
                    default=0.25,
                    type=float,
                    help="ewald factor(default: %(default)s)")

    parser.add_argument("-nx", "--nx", dest="nx",
                    default=8,
                    type=int,
                    help="nx grid points in x(default: %(default)s)")

    parser.add_argument("-ny", "--ny", dest="ny",
                    default=8,
                    type=int,
                    help="ny grid points in x(default: %(default)s)")

    parser.add_argument("-nz", "--nz", dest="nz",
                    default=18,
                    type=int,
                    help="nz grid points in x(default: %(default)s)")

    parser.add_argument('--create_job_flag', action='store_true', help='run actual analysis')

    parser.add_argument("-wm", "--wallminutes", dest="wallminutes",
                    default=50,
                    type=int,
                    help="minuts to run (default: %(default)s)")
    parser.add_argument("-wh", "--wallhours", dest="wallhours",
                    default=3,
                    type=int,
                    help="hours to run (default: %(default)s)")
    parser.add_argument("-pr", "--procs", dest="procs",
                    default=1,
                    type=int,
                    help="number of procs(default: %(default)s)")

    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    print(termcolor.colored('doing analysis', 'red'))
    run_analysis(args)
