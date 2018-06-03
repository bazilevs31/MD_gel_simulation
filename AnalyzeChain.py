#!/usr/bin/env python

# Program: AnalyzeChain.py
# Purpose: analyzes chain crystallinity
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeChain.py for help,
# example: python AnalyzeChain.py -mu 100 -s 20 -nch 100
# Requires: read_parameters.py, save_plots

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import MDAnalysis

import read_parameters
import save_plots
import AnalyzeLog
import get_path_names
import numba

def moving_average(a, n=3) :
    ret = np.cumsum(a, axis=0, dtype=float)
    ret[n:,:] = ret[n:,:] - ret[:-n,:]
    return ret[n - 1:,:] / n

# @numba.autojit
def get_chain_crystallinity(bonds, threshold, neigh):
    """
    input: bonds, neigh = number of neighbours to get the mean direction
    bonds = array of chain bonds
    gets chords of a chain
    then loops from within nleft:-nright of an array
    calculates the crystallinity of a chain

    crystallinity = aligned parts / all parts of a chain
    align part = cos2_local > 0.9
    """
    nleft = int(neigh)/2  # the left part of a chain
    nright = neigh - nleft
    # the right part of a chain which doesn't have any neighbours
    crystallinity = 0.
    cos2_local = 0.
    kk = 0.0
    Nbonds = len(bonds)
    AlignVec = moving_average(bonds, neigh)
    for i in range(nleft, Nbonds-nright):
        cos2_local = np.dot(AlignVec[i-nleft, :], bonds[i, :])
        cos2_local = 2.0*cos2_local**2-1.0
        if cos2_local > threshold:
            kk += 1.0  # if chain is aligned count the number of such chains
    if kk == 0:
        crystallinity = 0.0  # if there are no segments with aligned parts, i.e kk=0, then crystallinity == 0.
    else:
        crystallinity = kk/float(Nbonds-nleft-nright)
    return crystallinity

def get_bondlist_coords(u):
    """
    input: Universe
    output: bonds (that are in the domain, normalized)

    generate normalized coordinates of bond vectors
    get universe , return bonds(coordinates)
    generate coor of all bonds(bond = chord i-1 - i+1 ), normalize it
    """
    bonds = u.bonds.atom2.positions - u.bonds.atom1.positions
    # angles = u.angles
    # bonds = angles.atom3.positions - angles.atom1.positions
    # coords = angles.atom2.positions
    norm = np.linalg.norm(bonds, axis=1)
    bonds /= norm[:, None]  # the norm vector is a (nx1) and we have to
#   create dummy directions -> (n,3)
    return bonds


def calc_cryst(args):
    """
    Analyzes crystallinity
    """
    u = MDAnalysis.Universe(args.psffile, args.traj)
    psffile = get_path_names.get_filename(args.psffile)
    Nres = len(u.residues)
    frame = []
    g2 = []
    res_g2 = 0.
    # key parameter is the parameter used for parsing the files
    # keyparameter = get_path_names.get_filename(args.traj)
    keyparameter = args.keyparameter
    if args.filedumpskip:
        dumpskip = AnalyzeLog.get_dumpskip(filename='dumpskip.txt',
                                           keyparameter="dump" +
                                           keyparameter+"skip")
        timestep = AnalyzeLog.get_timestep(filename='dumpskip.txt',
                                           keyparameter="time" +
                                           keyparameter
                                           + "step")
    else:
        print "using given parameters"
        dumpskip = args.dumpskip
        timestep = args.timestep

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        all_g2 = 0.
        frame.append(ts.frame)
        for res in u.residues:
            chords = get_bondlist_coords(res)
            res_g2 = get_chain_crystallinity(chords,
                                             threshold=args.threshold,
                                             neigh=args.neigh)
            all_g2 += res_g2
        print ('frame {frame} , g2 = {cryst}'.format(frame=ts.frame,
                                                     cryst=all_g2/float(Nres)))
        g2.append(all_g2 / float(Nres))
    frame = np.array(frame)
    g2 = np.array(g2)
    save_plots.save_plot(frame, g2,
                         plotname='ICC',
                         name=psffile,
                         logplot=args.logplot,
                         dumpskip=dumpskip,
                         timestep=timestep,
                         duplicate=False)
    # save_plots.save_plot_cryst(frame, g2, psffile,
    #                            logplot=args.logplot,
    #                            dumpskip=dumpskip,
    #                            timestep=timestep)
    return None


def main():
    """main programm"""
    args = read_parameters.read_traj_vmd()
    calc_cryst(args)

if __name__ == '__main__':
    main()
