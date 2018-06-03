#!/usr/bin/env python


# Program: AnalyzeChain.py
# Purpose: analyzes counter ion condensation
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python Analyze_Condensation.py for help,

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import MDAnalysis as mda
import os

# import distances
import termcolor
import argparse
import pandas as pd
import Analyze_file_info



def find_boundary(u, f1, f2, box, args, verbose=False):
    """
    finds the boundary in the system

    Args:
        u (MDAnalysis Universe): MDanalysis trajectory
        f1 (float): ionization of the left part <= 1.
        f2 (float): ionization of the right part <=1.
        box (numpy array): box dimensions of the trajectory [lx, ly, lz]
        verbose (bool, optional): print the results

    Returns:
        right_min (float): z coordinate of the boundary
    """
    # f1 = 0.1
    # f2 = 1.
    # Nmon = 99

    left_max = 0
    right_min = box[2]

    if f2 == 0.3:
        f2 = 1/3.
    print(box)
    print('params {f1}, {f2}'.format(f1=f1, f2=f2))


    print('running condensation analysis')


    frames_array = []
    boundary_left_array = []
    boundary_right_array = []

    print("my traj parameters")
    print(args.startframe, args.endframe, args.dumpskip)
    # print("my trajectory with {nframes}".format(nframes=u.trajectory.n_frames))

    for ts in u.trajectory[args.startframe:args.endframe:args.dumpskip]:
    # for ts in u.trajectory[4:8:trajskip]:
        print('analyzing frame {0}'.format(ts.frame))
        frames_array.append(ts.frame)
        for r in u.residues[1:]:
            f = np.abs(r.charge/r.mass)
            mymax = r.atoms.positions[:,2].max()
            mymin = r.atoms.positions[:,2].min()
            if verbose:
                print('resid f = {f}, fr = {fr}, close? {close} min {mymin}, max {mymax} reslen {reslen}'.format(f=f,fr=f2, close=np.isclose(f,f2,0.05), mymin=mymin,mymax=mymax, reslen=r.atoms.n_atoms))
            if abs(mymax-mymin)<0.9*box[2]:
                if np.isclose(f,f2,0.05):
                    d = r.atoms.positions[:,2].min()
                    if d < right_min:
                        if verbose:
                            print('params {f2}, d right = {d}'.format(f2=f2, d=d))
                        right_min = d
                if np.isclose(f,f1,0.05):
                    d = r.atoms.positions[:,2].max()
                    if d > left_max:
                        if verbose:
                            print('params {f1}, dleft = {d}'.format(f1=f1, d=d))
                        left_max = d
        boundary_right_array.append(right_min)
        boundary_left_array.append(left_max)
    if verbose:
        print('params {f1}, {f2}, boundary theor = {theor}, left = {left}, f2 ioni boundary {right}'.format(f1=f1, f2=f2, theor=box[2]/2, left=left_max,right=right_min))
    return frames_array, boundary_right_array, boundary_left_array


def save_csv(temp_array, q_left_ave_array, q_right_ave_array):
    import pandas as pd
    df = pd.DataFrame()
    df['Temp'] = temp_array
    df['q_left'] = q_left_ave_array
    df['q_right'] = q_right_ave_array
    df['free_left'] = 1 - df['q_left']
    df['free_right'] = 1 - df['q_right']
    df.to_csv('countrion_cond.csv')



def analyze_ci_condensation(filename, tempname, args):
    """
    analyzes a dumy trajectory
    """
    analysis_name = '{filename}_t{tempname}'.format(filename=filename,tempname=tempname)
    datafile = '{filename}_t{tempname}{suffix}.data'.format(filename=filename,tempname=tempname, suffix=args.data_suffix)
    dcdfile = '{filename}_t{tempname}{suffix}.dcd'.format(filename=filename,tempname=tempname,suffix=args.traj_suffix)
    # u = mda.Universe('test_gel.data')
    f1, f2, _ , _ = Analyze_file_info.get_profile_params(analysis_name)
    # print(dcdfile)
    print("f1 = {f1}, f2={f2}".format(f1=f1, f2=f2))
    print('data',datafile)
    print('traj',dcdfile)
    print('it exists?')

    if (os.path.exists(datafile) and os.path.exists(dcdfile)):
        print('yes it does')
        print(os.path.exists(datafile))
        print(os.path.exists(dcdfile))

        u = mda.Universe(datafile,dcdfile)
        u.atoms.pack_into_box()
        box = u.trajectory.ts.dimensions[:3]
        # u = mda.Universe(datafile,dcdfile)

        print("analyzing traj with N_frames = {0}".format(u.trajectory.n_frames))

        frames_array, boundary_right_array, boundary_left_array = find_boundary(u,f1, f2, box, args)


        # # print("left = " + q_array_left)
        print("right")
        print(boundary_right_array)
        # print(q_array_left)
        print("left")
        print(boundary_left_array)


        # print(q_array_right)
        # ffP.clf()

        # return np.mean(q_array_left[2:]), np.mean(q_array_right[2:])
        # return np.mean(boundary_right_array[2:]), np.mean(boundary_left_array[2:])
        return frames_array, boundary_right_array, boundary_left_array
    else:
        print('doesnt exist')
        return None

def run_analysis(args):
    """runs the actual analysis"""

    files = args.files
    plot_flag = False
    # files = ['sc1844_nm100_l50_fl01_fr02','sc1844_nm100_l50_fl01_fr03','sc1844_nm100_l50_fl01_fr05','sc1844_nm100_l50_fl01_fr1','sc1844_nm100_l50_fl025_fr05']
    # file_mark_list = ['r','m','g','b','k']
    # file_marker_dict = dict(zip(files, file_mark_list))


    temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
    # for filename, _ in file_marker_dict.items():
    for filename in files:
        boundary_left_temp_array = []
        boundary_right_temp_array = []
        t_array = []

        if plot_flag:
            import matplotlib.pyplot as plt
            plt.style.use('seaborn')
            fig= plt.figure(1)
            fig.clf()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            fig2 = plt.figure(2)
            fig2.clf()
            axT = fig2.add_subplot(211)
            axTdelta = fig2.add_subplot(212)

        for tempname in temp_name_array:
            temp = Analyze_file_info.get_temp_from_string(tempname)
            # q1, q2 = 1, 2

            result_tuple = analyze_ci_condensation(filename, tempname, args)
            if result_tuple is not None:
                frames_array, boundary_right_array, boundary_left_array = result_tuple[0], result_tuple[1], result_tuple[2]
                if plot_flag:
                    ax1.plot(np.array(frames_array), np.array(boundary_right_array)/1000., '--',label='right T = {tempname}'.format(tempname=tempname))
                    ax2.plot(np.array(frames_array), np.array(boundary_left_array)/1000., '--',label='left T = {tempname}'.format(tempname=tempname))

                boundary_left_temp_array.append(np.mean(boundary_left_array[2:]))
                boundary_right_temp_array.append(np.mean(boundary_right_array[2:]))
                t_array.append(temp)
            # # i = df.loc[df['T'] == temp]
            # i = df.index[df['T'] == temp]
            # # mask = df['T'].values == temp
            # df['q1'][i] = q1
            # df['q2'][i] = q2
            # df.to_csv('q'+filename+'.csv')
        if plot_flag:
            ax1.set(xlim=[args.startframe, args.endframe],ylim=[0.2,0.65], ylabel='$z_{r}/L_z$')
            ax2.set(xlim=[args.startframe, args.endframe],ylim=[0.2,0.65], xlabel='frame', ylabel='$z_{l}/L_z$')
            axT.set(xlim=[0, 1.2],ylim=[0, 0.7],  ylabel='$z_{l}/L_z$')
            axTdelta.set(xlim=[0, 1.2],ylim=[-100, 400], xlabel='T, $[\\varepsilon]$', ylabel='$\Delta l$')

            axT.plot(np.array(t_array),np.array(boundary_right_temp_array)/1000.,'o--',label='right')
            axT.plot(np.array(t_array),np.array(boundary_left_temp_array)/1000.,'o--',label='left')
            ave_boundary_array = 0.5*(np.array(boundary_right_temp_array)+np.array(boundary_left_temp_array))
            error_boundary_array = np.abs(np.array(boundary_right_temp_array)-np.array(boundary_left_temp_array))
            delta_boundary_array = 9*50. - ave_boundary_array

            axTdelta.errorbar(np.array(t_array), delta_boundary_array,yerr= error_boundary_array,marker='o', c='green', ls='--')

            ax1.get_shared_x_axes().join(ax1, ax2)
            ax1.set_xticklabels([])

            axT.get_shared_x_axes().join(axT, axTdelta)
            axT.set_xticklabels([])

            ax1.legend()
            ax2.legend()
            axT.legend()

        df = pd.DataFrame([t_array, boundary_right_temp_array, boundary_left_temp_array])
        df = df.T
        df.columns = ['T', 'z1', 'z2']
        df['dl'] = 9*50 - 0.5*(df['z1'] + df['z2'])
        df.to_csv('{prefix}_{filename}.csv'.format(prefix=args.csv_prefix, filename=filename))
        if plot_flag:
            fig.savefig('{gelname}_boundary_analysis.png'.format(gelname=filename))
            fig2.savefig('{gelname}_boundary_analysis_temp.png'.format(gelname=filename))
    return None

# analyze_test_traj()
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
    parser.add_argument("-kf",
                    "--kframes",
                    dest="kframes",
                    default=9,
                    type=int,
                    help="Every k frames (ts.frame-1 / kframes  == 0) our cond_init array becomes np.ones(), which means that the counter ions have to stay within the backbone for kframes until they are considred condensed (default: %(default)s)")

    parser.add_argument("-r", "--rcut", dest="rcut",
                    default=2.,
                    type=float,
                    help="how close ci needs to be to the backbone (default: %(default)s)")

    parser.add_argument("-cs",
                        "--csv_prefix",
                        dest="csv_prefix",
                        default='Boundary',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")
    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    print(termcolor.colored('doing traj analysis', 'red'))
    run_analysis(args)

    # if args.test:
    #     print(termcolor.colored('doing traj test', 'red'))
    #     analyze_test_traj(args)
    # elif args.analysis:

    # elif args.testdata:
    #     print(termcolor.colored('doing test data', 'red'))
    #     analyze_test_data()

