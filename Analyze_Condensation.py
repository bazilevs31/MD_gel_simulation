#!/usr/bin/env python


# Program: AnalyzeChain.py
# Purpose: analyzes counter ion condensation
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python Analyze_Condensation.py for help,

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numba
import numpy as np
import MDAnalysis as mda
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import pandas as pd
# import distances
import termcolor
import argparse

import Analyze_file_info



# def create_slurm_string(procs, wh, dataname, tempname, temp):


@numba.autojit(nopython=True)
def calc_condensation(Gel_pos, CI_pos, CI_cond_array,box, rcut):
    """
    loops through each monomer of the gel
    looks for neighboring counter ions
    finds the counterions that are not yet flagged condensed CI_cond_array[i]=0
    and if they are close enough < rcut marks them as condensed
    CI_cond_array[i] = 1

    Args:
        Gel_pos (numpy array): positions of the gel monomers
        CI_pos (numpy array): positions of the counter ions
        CI_cond_array (numpy array): condensation index array per each ci
        rcut (float): distance at which the counter ion is considered condensed usually < 2.5

    Returns:
        CI_cond_array: Updated condensation index array per each ci
    """
    nci = len(CI_pos)
    ngel = len(Gel_pos)
    d2 = 0.0
    for ig in range(ngel):
        # print('gel pos ',g)
        for ici in range(nci):
            # print('ci pos ',ci)
            if CI_cond_array[ici] == 0:
                # d = dist_w_pbc(CI_pos[ici], Gel_pos[ig],box)
                X  = CI_pos[ici, 0] - Gel_pos[ig,0]
                Y  = CI_pos[ici, 1] - Gel_pos[ig,1]
                Z  = CI_pos[ici, 2] - Gel_pos[ig,2]

                # Periodic boundary condition
                X  -= box[0] * np.rint(X/box[0])
                Y  -= box[1] * np.rint(Y/box[1])
                Z  -= box[2] * np.rint(Z/box[2])
                d2 = X*X + Y*Y + Z*Z
                if d2 < rcut*rcut:
                    CI_cond_array[ici] += 1
    return CI_cond_array


def find_boundary(u, f1, f2, box, verbose=True):
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

    for r in u.residues:
        f = np.abs(r.charge/r.mass)
        mymax = r.atoms.positions[:,2].max()
        mymin = r.atoms.positions[:,2].min()
        if verbose:
            print('resid f = {f}, fr = {fr}, close? {close} min {mymin}, max {mymax} reslen {reslen}'.format(f=f,fr=f2, close=np.isclose(f,f2,0.05), mymin=mymin,mymax=mymax, reslen=r.atoms.n_atoms))
        if abs(mymax-mymin)<0.8*box[2]:
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
    if verbose:
        print('params {f1}, {f2}, boundary theor = {theor}, left = {left}, f2 ioni boundary {right}'.format(f1=f1, f2=f2, theor=box[2]/2, left=left_max,right=right_min))
    return right_min

def find_condensed_array(u, f1, f2, args):
    """analyzes the universe to find the number of condensed
    on the right and left

    Args:
        u (MDAnalysis Universe): the universe object with the trajectory
        f1
        f2
        args
        # k_frames=0, startframe=2, endframe=40, trajskip=4, rcut=2.
    Returns:
        q_array_left: numpy array with the fraction of condensed counterions
        on the left side of the box
        q_array_right: numpy array with the fraction of condensed counterions
        on the right side of the box
        frames_array: numpy of frames
    """
    # k_frames=0, startframe=2, endframe=40, trajskip=4, rcut=2.
    k_frames=args.kframes
    startframe=args.startframe
    endframe=args.endframe
    trajskip=args.dumpskip
    rcut=args.rcut

    q_left_array = []
    q_right_array = []

    # u.atoms.pack_into_box()

    box = u.trajectory.ts.dimensions[:-3]
    print("mybox")
    print(box)

    if args.analysis:
        print("analyzing boundary")
        boundary_pos = find_boundary(u, f1=f1, f2=f2, box=box)
    else:
        # boundary_pos = box[2]/2.
        pbc_side_pad = 100
        boundary_pos = (box[2] - pbc_side_pad)/2. - 10.
        # boundary_pos = 0.8*box[2]/2.


    gel_right = u.select_atoms('type 1 2 3 and prop z > {boundary}'.format(boundary=boundary_pos))
    cions_right = u.select_atoms('type 4 and prop z > {boundary}'.format(boundary=boundary_pos))

    print("right cions {0}".format(cions_right.n_atoms))

    gel_left = u.select_atoms('type 1 2 3 and prop z <= {boundary}'.format(boundary=boundary_pos))
    cions_left = u.select_atoms('type 4 and prop z <= {boundary}'.format(boundary=boundary_pos))

    print("left cions {0}".format(cions_left.n_atoms))




    gel_right_positions = gel_right.atoms.positions
    cions_right_positions = cions_right.atoms.positions
    gel_right_indices = gel_right.indices
    cions_right_indices = cions_right.indices

    gel_left_positions = gel_left.atoms.positions
    cions_left_positions = cions_left.atoms.positions
    gel_left_indices = gel_left.indices
    cions_left_indices = cions_left.indices

    Nci_right = cions_right.n_atoms
    Nci_left = cions_left.n_atoms

    ci_cond_array_right = np.zeros(Nci_right)
    ci_cond_array_left = np.zeros(Nci_left)

    print('running condensation analysis')


    frames_array = []
    q_array_left = []
    q_array_right = []

#
    ci_cond_array_right_init = np.ones(cions_right.n_atoms)
    ci_cond_array_left_init = np.ones(cions_left.n_atoms)
    print("my traj parameters")
    print(startframe, endframe, trajskip, k_frames)
    # print("my trajectory with {nframes}".format(nframes=u.trajectory.n_frames))

    for ts in u.trajectory[startframe:endframe:trajskip]:
    # for ts in u.trajectory[4:8:trajskip]:
        print('analyzing frame {0}'.format(ts.frame))
        frames_array.append(ts.frame)

        gel_right_positions = np.copy(gel_right.atoms.positions)
        cions_right_positions = np.copy(cions_right.atoms.positions)
        # gel_right_indices = gel_right.indices
        # cions_right_indices = cions_right.indices

        gel_left_positions = np.copy(gel_left.atoms.positions)
        cions_left_positions = np.copy(cions_left.atoms.positions)

        # print("gel left pos")
        # print(gel_left_positions)
        # gel_left_indices = np.copy(gel_left.indices)
        # cions_left_indices = np.copy(cions_left.indices)

        ci_cond_array_right = np.zeros(cions_right.n_atoms)
        ci_cond_array_left = np.zeros(cions_left.n_atoms)


        if (ts.frame-1) % (k_frames) == 0:
            print('tsframe-1%k_frames == 0')
            print(ts.frame)
            print(k_frames)
            q_index = np.count_nonzero(ci_cond_array_right_init>=1.)
            q_array_right.append(q_index)
            print("q right index = {q}".format(q=q_index/Nci_right))

            q_index = np.count_nonzero(ci_cond_array_left_init>=1.)
            q_array_left.append(q_index)

            print("q left index = {q}".format(q=q_index/Nci_left))

            ci_cond_array_right_init = np.ones(cions_right.n_atoms)
            ci_cond_array_left_init = np.ones(cions_left.n_atoms)

        calc_condensation(gel_right_positions, cions_right_positions, ci_cond_array_right,box, rcut)
        calc_condensation(gel_left_positions, cions_left_positions, ci_cond_array_left,box, rcut)
        # print('current condensation array')
        # print(np.unique(ci_cond_array_right))
        # print(np.unique(ci_cond_array_left))
        ci_cond_array_right_init *= ci_cond_array_right
        ci_cond_array_left_init *= ci_cond_array_left
        # print(ci_cond_array_right)



        # print("q index = {q}".format(q=q_index))

    # print(frames_array, q_array_right)
    # print(frames_array, q_array_left)

    q_array_right = np.array(q_array_right)/Nci_right
    q_array_left = np.array(q_array_left)/Nci_left

    # print("right")
    # print(q_array_right)

    # print("left")
    # print(q_array_left)

    return q_array_right,q_array_left, frames_array

def save_csv(temp_array, q_left_ave_array, q_right_ave_array):
    import pandas as pd
    df = pd.DataFrame()
    df['Temp'] = temp_array
    df['q_left'] = q_left_ave_array
    df['q_right'] = q_right_ave_array
    df['free_left'] = 1 - df['q_left']
    df['free_right'] = 1 - df['q_right']
    df.to_csv('countrion_cond.csv')


def analyze_test_data(N_repeat_array=[1, 4, 7], N_every_1_array=[3, 5,9, 15], N_every_2_array=[3, 5,9, 15]):

    df = pd.DataFrame()

    q_left_ave_array = []
    q_right_ave_array = []
    temp_name_array = ['01','03','11']
    temp_array = [0.1, 0.3, 1.1]
    # for rcut in [2.0, 3.5, 5.0, 10.]:
    for N_repeat in N_repeat_array:
        for N_every_1 in N_every_1_array:
            for N_every_2 in N_every_2_array:
                # if N_every_2 < N_every_1:
                #     continue
                datafile = "test_gel_long_{nevery1}_{nevery2}_r{nrepeat}.data".format(nevery1=N_every_1,nevery2=N_every_2, nrepeat=N_repeat)
                # tempname = '01'
                filename = 'test_gel_long_{nevery1}_{nevery2}_r{nrepeat}'.format(nevery1=N_every_1,nevery2=N_every_2, nrepeat=N_repeat)

                u = mda.Universe(datafile)
                print(filename)
                q_array_right, q_array_left, frames_array = find_condensed_array(u)

                # print("left = " + q_array_left)
                print("left")
                print(q_array_left)
                print("right")
                print(q_array_right)

def analyze_test_traj(args):
    """
    analyzes a dumy trajectory
    """
    datafile = "./test/new.data"
    dcdfile = "./test/new.dcd"
    u = mda.Universe(datafile,dcdfile)
    print(u.trajectory.n_frames)
    # print(u.trajectory)

    q_array_right, q_array_left, frames_array = find_condensed_array(u,f1=0.1, f2=1.0, args=args)

    # print("left = " + q_array_left)
    print("left")
    print(q_array_left)
    print("right")
    print(q_array_right)
    return None

def analyze_ci_condensation(filename, tempname, args):
    """
    analyzes a dumy trajectory
    """
    analysis_name = '{filename}_t{tempname}'.format(filename=filename,tempname=tempname)
    datafile = '{filename}_t{tempname}{suffix}.data'.format(filename=filename,tempname=tempname, suffix=args.data_suffix)
    dcdfile = '{filename}_t{tempname}{suffix}.dcd'.format(filename=filename,tempname=tempname,suffix=args.traj_suffix)f
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

        # u = mda.Universe(datafile,dcdfile)

        print("analyzing traj with N_frames = {0}".format(u.trajectory.n_frames))

        q_array_right, q_array_left, frames_array = find_condensed_array(u,f1, f2, args)

        # print("left = " + q_array_left)
        print("left")
        print(q_array_left)
        print("right")
        print(q_array_right)

        return np.mean(q_array_left[2:]), np.mean(q_array_right[2:])
    else:
        print('doesnt exist')
        return np.NaN, np.NaN

def run_analysis(args):
    """runs the actual analysis"""

    files = args.files
    # files = ['sc1844_nm100_l50_fl01_fr02','sc1844_nm100_l50_fl01_fr03','sc1844_nm100_l50_fl01_fr05','sc1844_nm100_l50_fl01_fr1','sc1844_nm100_l50_fl025_fr05']
    # file_mark_list = ['r','m','g','b','k']
    # file_marker_dict = dict(zip(files, file_mark_list))
    temp_name_array = ['01','02','04','04','05','06','07','08','09','10','11']

    # for filename, _ in file_marker_dict.items():
    for filename in files:
        q1_array = []
        q2_array = []
        t_array = []
        df = pd.DataFrame(columns = ['T', 'q1', 'q2'])
        # df.columns = ['T', 'q1', 'q2']
        df['T'] = np.linspace(0.1, 1.1, 11)

        for tempname in temp_name_array:
            temp = Analyze_file_info.get_temp_from_string(tempname)
            # q1, q2 = 1, 2
            q1, q2 = analyze_ci_condensation(filename, tempname, args)
            q1_array.append(q1)
            q2_array.append(q2)
            t_array.append(temp)
            # # i = df.loc[df['T'] == temp]
            # i = df.index[df['T'] == temp]
            # # mask = df['T'].values == temp
            # df['q1'][i] = q1
            # df['q2'][i] = q2
        df = pd.DataFrame([t_array, q1_array, q2_array])
        df = df.T
        df.columns = ['T', 'q1', 'q2']
        df.to_csv('{prefix}_{filename}.csv'.format(prefix=args.csv_prefix, filename=filename))
        # df.to_csv('q'+filename+'.csv')
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

    group = parser.add_mutually_exclusive_group()
    # group.add_argument('--lines', action='store_true', help='Plot the data with the lines style')
    group.add_argument('--test', action='store_true', help='run test analysis')
    group.add_argument('--testdata', action='store_true', help='run test analysis')
    group.add_argument('--analysis', action='store_true', help='run actual analysis')

    opts, rem_args = parser.parse_known_args()
    if opts.test:
        # use options and namespace from first parsing
        parser.add_argument("--testdatafile", dest="testdatafile",
                            default="./test/new.data",
                            type=lambda x: is_valid_file(parser, x),
                            help="Input trajectory file)", metavar="FILE")
        parser.add_argument("--testdcdfile", dest="testdcdfile",
                            default="./test/new.dcd",
                            type=lambda x: is_valid_file(parser, x),
                            help="Input trajectory file)", metavar="FILE")


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
                        default='COND_v1',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")
    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))


    if args.test:
        print(termcolor.colored('doing traj test', 'red'))
        analyze_test_traj(args)
    elif args.analysis:
        print(termcolor.colored('doing traj analysis', 'red'))
        run_analysis(args)
    elif args.testdata:
        print(termcolor.colored('doing test data', 'red'))
        analyze_test_data()

