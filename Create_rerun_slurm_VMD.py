#!/usr/bin/env python

import termcolor
import argparse
import os
import shutil

def create_run_sh(**params):
    """creates the run sh
    with sbatch name.sh
    """
    with open('run.sh', 'a') as f:
        f.write("\n sbatch {0[slurm_out_name]}.sh \n".format(params))
    return None

# def create_slurm_string(procs, wh, dataname, tempname, temp):
def create_slurm_header(**params):
    """
    creates a slurm string header
    that will be later populated by the lammps or vmd rerun/analyze
    commands
    """
    flag_vmd = params['flag_vmd']
    flag_rerun_lammps = params['flag_rerun_lammps']

    if flag_vmd:
        params['mem_string'] = '--mem=82G'
        module_string = '\n ##writing a VMD dx analyze file \n'
        module_string += '\n module load vmd/1.9.3 \n '

    elif flag_rerun_lammps:
        params['mem_string'] = '--mem-per-cpu 600'
        module_string = '\n ##writing a Lammps rerun file \n'
        module_string += '\n module load intel/2016.4 openmpi/2.1.1 \n '
        module_string += '\n module load nixpkgs/16.09 lammps-user-intel/20170811 \n '
    else:
        raise RuntimeError('there should be either flag_vmd or flag_rerun_lammps')

    slurm_template_file = 'template_slurm.sh'
    with open(slurm_template_file) as f:
        read_data = f.read()
        heading_script = (read_data.format(params))

    with open("{0[slurm_out_name]}.sh".format(params), 'w') as f:
        f.write(heading_script  + module_string +  '\n' +'\n')

    return heading_script + module_string

def add_slurm_run_commands(**params):

    flag_vmd = params['flag_vmd']
    flag_rerun_lammps = params['flag_rerun_lammps']


    if flag_vmd:
        runscript = 'vmd -dispdev text -e {0[vmd_out_name]}.vmd'.format(params)

    elif flag_rerun_lammps:
        runscript = '\n srun lmp_intel_cpu_openmpi -in {0[rerun_name]}.lammps -echo both \n'.format(params)
    else:
        raise RuntimeError('there should be either flag_vmd or flag_rerun_lammps')


    with open("{0[slurm_out_name]}.sh".format(params), 'a') as f:
        f.write(runscript + '\n' +'\n')

    return runscript

def add_vmd_run_commands(**params):
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
    pmepot -mol {0[molnumber]} -ewaldfactor {0[ewaldfactor]}  -frames {0[startframe]}:{0[trajskip]}:{0[endframe]}  -grid {{ {0[nx]} {0[ny]} {0[nz]} }}    -dxfile {0[out_folder]}/{0[dxname]}.dx -loadmol none
    mol delete {0[molnumber]}

    """.format(params)


    with open("{0[vmd_out_name]}.vmd".format(params), 'a') as f:
        f.write(headstring  + body_string +  '\n' +'\n')
    return headstring  + body_string

def create_lammps_rerun_file(**params):
    """
    creates a unique lammps rerun file
    for a certain temperature and filenam
    """
    # params_dict = vars(args)

    # params = {
    # 'filename' : 'init',
    # 'out_folder' : 'profiles_v2',
    # 'datafile' : 'init.data',
    # 'data_folder' : 'data',
    # 'dcdfile' : 'init.dcd',
    # 'dcd_folder' : 'data',
    # 'Dx' : 0.11,
    # 'trajectory_skip' : 10,
    # 'lammps_out_name' : 'rerun_profiles'
    # }

    rerun_template_file_init = 'template_rerun_init.lammps'
    rerun_template_file_profiles = 'template_rerun_profiles.lammps'

    with open(rerun_template_file_init) as f:
        read_data = f.read()
        init_string = (read_data.format(params))

    with open(rerun_template_file_profiles) as f:
        read_data = f.read()
        profiles_string = (read_data.format(params))

    result_string = '\n' + init_string +'\n'+ profiles_string + '\n'
    print(result_string)

    with open('{0[rerun_name]}.lammps'.format(params), 'w') as f:
        f.write(result_string)

    return result_string

def run_analysis(args):
    """runs the actual analysis"""

    params_dict = vars(args)

    print('==============')
    print('please copy all the template files in order for this program to work well, for example you could use')
    print('cp ~/packages/Python_packages/template_* .')
    print('==============')

    # shutil.copyfile('/home/bazilevs/packages/Python_packages/template*', './')

    with open("{0[vmd_out_name]}.vmd".format(params_dict), 'w') as f:
        f.write('\n')

    temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
    molnumber = 0
    # for filename, _ in file_marker_dict.items():
    create_slurm_header(**params_dict)

    print("files")
    print(params_dict['files'])
    for filename in params_dict['files']:
        print("analyzing file")
        print(filename)
        for tempname in temp_name_array:
            params_dict['filename'] = filename
            params_dict['tempname'] = tempname
            dcdfile = '{0[dcd_folder]}/{0[filename]}_t{0[tempname]}_{0[dcd_suffix]}.dcd'.format(params_dict)
            datafile = '{0[data_folder]}/{0[filename]}_t{0[tempname]}_{0[data_suffix]}.data'.format(params_dict)
            dxname = '{0[filename]}_t{0[tempname]}'.format(params_dict)
            rerun_name = '{0[filename]}_t{0[tempname]}'.format(params_dict)

            if os.path.isfile(dcdfile) and os.path.isfile(datafile):
                print('file exists')
                params_dict['datafile'] = datafile
                params_dict['dcdfile'] = dcdfile
                params_dict['dxname'] = dxname
                params_dict['molnumber'] = molnumber
                params_dict['rerun_name'] = rerun_name
                print(params_dict)
                print(termcolor.colored('adding run commands with the parameters', 'blue'))
                print(params_dict)
                if params_dict['flag_rerun_lammps']:
                    create_lammps_rerun_file(**params_dict)
                    add_slurm_run_commands(**params_dict)
                    create_run_sh(**params_dict)
                elif params_dict['flag_vmd']:
                    add_vmd_run_commands(**params_dict)
                # write_lammps_to_slurm(**params_dict)
                molnumber += 1

        # finishing the writing with the exit
    if params_dict['flag_vmd']:
        add_slurm_run_commands(**params_dict)

    with open("{0[vmd_out_name]}.vmd".format(params_dict), 'a') as f:
        f.write(' \n exit \n')

    # if params_dict['create_job_flag']:
    #     with open("slurm_{0[vmd_out_name]}.sh".format(params_dict), 'w') as f:
    #         f.write(create_slurm_string(**params_dict))
    #         f.write(' \n exit \n')

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

    group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('--lines', action='store_true', help='Plot the data with the lines style')
    group.add_argument('--flag_vmd', action='store_true', help='run test analysis')
    group.add_argument('--flag_rerun_lammps', action='store_true', help='run test analysis')

    parser.add_argument('--lammps_potential_name',dest="lammps_potential_name", default='total', help='what potential do you want? could be total-coul-lj')


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


    parser.add_argument("-ol",
                        "--slurm_out_name",
                        dest="slurm_out_name",
                        default='rerun_profiles',
                        type=str,
                        help="the name of output lamps (default: %(default)s)")


    parser.add_argument("-Dx", "--Dx", dest="Dx",
                    default=0.05,
                    type=float,
                    help="the size of the bin if the total length is 1. (default: %(default)s)")


    parser.add_argument("-e", "--endframe", dest="endframe",
                    default=30,
                    type=int,
                    help="End frame of the trajectory file type (default: %(default)s)")
    parser.add_argument("-st",
                    "--startframe",
                    dest="startframe",
                    default=2,
                    type=int,
                    help="Start frame of the trajectory file type (default: %(default)s)")

    parser.add_argument("-ts",
                        "--trajskip",
                        dest="trajskip",
                        default=3,
                        type=int,
                        help="This is a factor for transforming trajectory frames into lj units, numframes*dumpskip = simtime [lj] (default: %(default)s)")
    # parser.add_argument('--create_job_flag', action='store_true', help='run actual analysis')



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

    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    print(termcolor.colored('doing analysis', 'red'))
    run_analysis(args)
