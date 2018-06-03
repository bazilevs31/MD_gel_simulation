#!/usr/bin/env python

import pandas as pd
from gridData import Grid
import glob
import os.path
import Analyze_file_info
import argparse
import termcolor

def find_dx_file(dirname, filename):
    print( "looking for matches for the file", filename)

    _, gelname, _  = Analyze_file_info.get_profile_filename(filename)
    # gelname = filename
    for fname in glob.glob(dirname+'/*.dx'):
        print( "looking for dx files")
        # analyze_Pressure_profile("energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt")
        if os.path.isfile(fname) and (gelname in fname):
            print( "file exists", "i found it gelnmae, fname, filename ",gelname, fname, filename)
            return fname, gelname
    return None

def create_csv_files(**params):
    """goes through the directory and finds files matching the press param and then analyzes the pressure and dx profiles, returns a csv file where all of this information is neatly complied in a coherent fashion"""

    # dirname = '/Users/bazilevs/Work/cedar/18-02-05'
    # for fname in glob.glob(dirname+'/profiles/pressalltot*.txt'):
    # dirname  = args.dirname
    # dx_dirname = args.dx_folder
    # press_search_param = args.profilename
    # press_dirname = args.profiles_folder

    for fname in glob.glob('{0[profiles_folder]}/{0[profilename]}*.txt'.format(params)):

    # for fname in glob.glob(dirname+'/profiles/pressalltot*.txt'):
        print( fname)
        # analyze_Pressure_profile("energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt")
        # dx_ = '{dirname}/{dx_dirname}'.format(dirname=dirname, dx_dirname=dx_dirname)
        # _tmp = os.path.basename(fname)
        # _tmp_filename = _tmp.replace('{press_search_param}_'.format(press_search_param=press_search_param),'')
        # _out = find_dx_file(dxfile_with_path, _tmp_filename.replace('.txt',''))
        _out = find_dx_file(params['dx_folder'], os.path.basename(fname))
        if _out is not None:
            dxfile, gelname = _out

        if os.path.isfile(fname) and (_out is not None):
            print( "file exists now analyzing the pressure")
            print( fname)
            print( "looking for dx files")

            df_p = Analyze_file_info.analyze_pressure_profile(fname)
            df_phi = Analyze_file_info.analyze_dx(dxfile)

            df_means = df_p.groupby(df_p['Id']).mean()
            df_stds = df_p.groupby(df_p['Id']).std()
            df_stds.add_suffix('_Error')
            df_p_full = pd.concat([df_means, df_stds], axis=1)

            print( df_p.head())
            print( df_phi.head())
            # result = pd.concat([df_p, df_means, df_stds], axis=1)

            # ci_name =
            ci_name = fname.replace('alltot','cipair')
            # ci_name = _.replace('total','coul')
            # ci_name = fname

            df_ci = Analyze_file_info.analyze_pressure_profile(ci_name)
            df_ci_means = df_ci.groupby(df_ci['Id']).mean()
            df_ci_stds = df_ci.groupby(df_ci['Id']).std()
            print("cion file is {ci_name}".format(ci_name=ci_name))
            print("====== CIONS ======")
            print(df_ci_means)

            result = pd.concat([df_phi, df_means,df_stds.add_suffix('_error'),df_ci_means.add_suffix('_ci'),df_ci_stds.add_suffix('_ci_error')], axis=1)

            result.to_csv('./{gelname}.csv'.format(gelname=gelname))
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

# dirname = '/Users/bazilevs/Work/cedar/18-03-19'
# dirname = os.path.abspath(dirname)
# press_dirname = 'profiles'
# press_search_param = 'pressalltot'
# dx_dirname = 'dx'

    # parser.add_argument('-d', "--dirname", dest="dirname",
    #                         default='/Users/bazilevs/Work/cedar/18-03-19',
    #                         type=lambda x: is_valid_file(parser, x),
    #                         help="folder where profiles are)", metavar="FILE")

    # parser.add_argument('-xdir', "--dx_dirname", dest="dx_dirname",
    #                         default='/Users/bazilevs/Work/cedar/18-03-19',
    #                         type=lambda x: is_valid_file(parser, x),
    #                         help="folder where dx folder is are)", metavar="FILE")

    parser.add_argument("-dx",
                        "--dx_folder",
                        dest="dx_folder",
                        default='./dx',
                        type=str,
                        help="default folder for dx (default: %(default)s)")

    parser.add_argument("-pd",
                        "--profiles_folder",
                        dest="profiles_folder",
                        default='./profiles',
                        type=str,
                        help="default folder for dx (default: %(default)s)")

    parser.add_argument("-pn",
                        "--profilename",
                        dest="profilename",
                        default='pressalltot',
                        type=str,
                        help="the name of profilename (default: %(default)s)")
    parser.add_argument("-pci",
                        "--ci_profile_name",
                        dest="ci_profile_name",
                        default='presscitot',
                        type=str,
                        help="the name of profilename (default: %(default)s)")

    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    print(termcolor.colored('doing analysis', 'red'))
    create_csv_files(**vars(args))
