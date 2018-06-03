#!/usr/bin/env python


import MDAnalysis as mda
import MDAnalysis.analysis.lineardensity
import os

# import distances
import termcolor
import argparse


def analyze_traj(args, filename, tempname):
    """analyzes the lienar density for a given universe and save the output files"""
    # u =
    datafile = '{filename}_t{tempname}{suffix}.data'.format(filename=filename,tempname=tempname, suffix=args.data_suffix)
    dcdfile = '{filename}_t{tempname}{suffix}.dcd'.format(filename=filename,tempname=tempname,suffix=args.traj_suffix)
    # u = mda.Universe('test_gel.data')
    print('data',datafile)
    print('traj',dcdfile)
    print('it exists?')

    if (os.path.exists(datafile) and os.path.exists(dcdfile)):
        print('yes it does')
        print(os.path.exists(datafile))
        print(os.path.exists(dcdfile))

        u = mda.Universe(datafile,dcdfile)
        u.atoms.pack_into_box()
        # u = mda.Universe(datafile,dcdfile)

        print("analyzing traj with N_frames = {0}".format(u.trajectory.n_frames))
        Lx, Ly, Lz = u.trajectory.ts.dimensions[:3]
        Nbins = args.nbins
        dz = Lz/float(Nbins)

        gel = u.select_atoms('type 1 2 3')
        cions = u.select_atoms('type 4')
        # condensed = cions.select_atoms('around 3.5 global group polymer', polymer=gel, updating=True)
        ldens = MDAnalysis.analysis.lineardensity.LinearDensity(gel,grouping='atoms', binsize=dz,start=2)
        ldens.run()
        ldens.save(description='gel_n{nbins}'.format(nbins=args.nbins), form='txt')

        ldens = MDAnalysis.analysis.lineardensity.LinearDensity(cions,grouping='atoms', binsize=dz,start=2)
        ldens.run()
        ldens.save(description='cions_n{nbins}'.format(nbins=args.nbins), form='txt')
    else:
        print('doesnt exist')
    return None


def run_analysis(args):
    """runs analysis with the given args"""
    # u = mda.Universe('sc1844_nm100_l50_fl01_fr1_t01_ljeq1_coul_eq_coulsim.data',
        # 'sc1844_nm100_l50_fl01_fr1_t01_ljeq1_coul_eq_coulsim_press.dcd')

    temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
    # for filename, _ in file_marker_dict.items():
    for filename in args.files:
        for tempname in temp_name_array:
            analyze_traj(args, filename, tempname)
    return None
# selection, , **kwargs

if __name__ == '__main__':
    # main()


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


    parser.add_argument('--analysis', action='store_true', help='run actual analysis')

    parser.add_argument('-a', action='append', dest='files',
                        type=str,
                        default=[],
                        help='Add values to analyze values to a list type (default: %(default)s)',
                        )
    parser.add_argument("-datas",
                        "--data_suffix",
                        dest="data_suffix",
                        default='_ljeq1_coul_eq_coulsim',
                        type=str,
                        help="suffix for the datafile gelname_<suffix> (default: %(default)s)")
    parser.add_argument("-trs",
                        "--traj_suffix",
                        dest="traj_suffix",
                        default='_ljeq1_coul_eq_coulsim_press',
                        type=str,
                        help="suffix for the dcdfile gelname_<suffix> (default: %(default)s)")

    # k_frames=9,startframe=2, endframe=30,trajskip=3

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

    parser.add_argument("-ds",
                        "--dumpskip",
                        dest="dumpskip",
                        default=3,
                        type=int,
                        help="This is a factor for transforming trajectory frames into lj units, numframes*dumpskip = simtime [lj] (default: %(default)s)")

    parser.add_argument("-n", "--nbins", dest="nbins",
                    default=15,
                    type=int,
                    help="how many beans do we need (default: %(default)s)")

    parser.add_argument("-p",
                        "--postfix",
                        dest="postfix",
                        default='nomidboundary',
                        type=str,
                        help="prefix for the output  (default: %(default)s)")
    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    print(termcolor.colored('doing traj analysis', 'red'))
    run_analysis(args)
