#!/usr/bin/env python

# Program: AnalyzeLamellae.py
# Purpose: calculates number of lamellae/all chains in the system
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeLamellae.py -h for help,
# Requires: read_parameters.py, save_plot

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


# todo: maybe we should normalize the number of lamellaes by number of crystallized chains, rather
# than number of all bonds
import os
import numpy as np
# from MDAnalysis import *
import MDAnalysis
# from read_parameters import read_traj_vmd
import read_parameters
import save_plots
import termcolor
import AnalyzeChain
import get_path_names
import AnalyzeLog
import numba
# this alghoirthm looks for the lamellae detection:

# @numba.jit
def get_seg_cryst(i, bonds, threshold=0.7, neigh=4):
    nleft, nright = int(neigh)/2, int(neigh)/2
    cos2_local = 0.
    AlignVec = AnalyzeChain.moving_average(bonds, neigh)
    cos2_local = np.dot(AlignVec[i-nleft, :], bonds[i, :])
    # cos2_local += np.dot(AlignVec[i+nright, :], bonds[i, :])
    # cos2_local /= 2.
    P2 = (3.*cos2_local**2.0 - 1.)/2.
    return P2


def main(args):
    """main func"""
    psffile = get_path_names.get_filename(args.psffile)
    u = MDAnalysis.Universe(psffile+'.psf', args.traj)

    bonds = AnalyzeChain.get_bondlist_coords(u.residues[0])
    L = 0.  # number of lamellae
    # the skip on the left and on the right,
    nleft, nright = int(args.neigh)/2, int(args.neigh)/2
    Nbonds = len(bonds)
    Nres = len(u.residues)
    timeframes, L_array = [], []

    if args.filedumpskip:
        dumpskip = AnalyzeLog.get_dumpskip(filename='dumpskip.txt',
                                           keyparameter="dump" +
                                           args.keyparameter+"skip")
        timestep = AnalyzeLog.get_timestep(filename='dumpskip.txt',
                                           keyparameter="time" +
                                           args.keyparameter
                                           + "step")
    else:
        print "using given parameters"
        dumpskip = args.dumpskip
        timestep = args.timestep

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        print "frame = {frame}, percent of lamellae={num_of_lam}"\
               .format(frame=ts.frame,
                       num_of_lam=L/float(Nres))
        timeframes.append(ts.frame)
        L = 0.
        for myres in u.residues:
            bonds = AnalyzeChain.get_bondlist_coords(myres)
            Flag = True  # assume we hit the straight segment of the chain
            atom_indices = []  # indices of atoms in the current lamellae
            for i in xrange(nleft, Nbonds-nright-1):
                if not Flag:  # if we ended a straight segment, check the len
                    Nlam = len(atom_indices)
                    if Nlam > 9:
                        L += 1.
# print colored( 'lamellae size is %d ' % Nlam, 'green')
                    # start looking for a new lamellae
                    atom_indices = []
                    Flag = True
                else:
                    # if we still are hitting straight segments
                    ci = get_seg_cryst(i, bonds,
                                       threshold=args.threshold,
                                       neigh=args.neigh)
                    if ci > args.threshold:
                        # print "segment = %d, | ci = %4.3f" % (i,ci)
                        Flag = True
                        atom_indices.append(i)
                    elif ci > 0.0:
                        # print "segment = %d, > ci = %4.3f" % (i,ci)
                        Flag = False
                    else:
                        ValueError("crystallinity is negative")
        L_array.append(L/float(Nres))
    timeframes = np.array(timeframes)
    L_array = np.array(L_array)
    save_plots.save_plot(timeframes, L_array,
                         plotname='lamellae',
                         name=psffile,
                         logplot=args.logplot,
                         dumpskip=dumpskip,
                         timestep=timestep,
                         ylabel='frac{N_{stems}}{N_{chains}}')

    # print "number of lammelae in the chain 1 is = %f" % float(L/Nres)

if __name__ == '__main__':

    args = read_parameters.read_traj_vmd()
    main(args)
