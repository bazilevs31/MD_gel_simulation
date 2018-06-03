#!/usr/bin/env python

# Program: AnalyzeEnd2End.py
# Purpose: Calculates end-to-end distance, can handle polydisperse systems
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeEnd2End.py -h for help, 
# Requires: read_parameters.py, save_plots.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 

import numpy as np
from MDAnalysis import *
import os
from read_parameters import read_traj_vmd
from save_plots import save_plot_r2n

def max_chain(u):
    """ 
    input: MDAnalysis universe 
    output: maximum length of all of the chains present(integer)
    """
    maxlen=0
    for res in u.residues:
        reslen=len(res)
        if maxlen<reslen:
            maxlen=reslen
    # print "maxlen = %f" % maxlen
    return  maxlen-1


def get_r2n(u,psffile,frame,Noff=1,logplot=False):
    """
    create a list of 
    R2_array - array of distances
    n_array - array of number of bonds between atoms
    k_array - array of number of atoms with this bonds
    start looping in residues
    for every residue:
        start from the middle of it 
            calculate the closest atom_i - atom_-i
            if it is the first time we have the number of bonds so big, we expand our lists by appending
            else: we just put it to the nth position
    then the last elements of the array will be deleted, since there is not enough statistics for this chains

    """
    R2_array = []
    n_array = []
    k_array = []

    for res in u.residues:
        chainlen = len(res)/2
        for i in range(chainlen)[::-1]:
            ag1, ag2 = res.atoms[i].pos, res.atoms[-i-1].pos
            tmpdiff = ag1-ag2
            r2 = np.dot(tmpdiff,tmpdiff)
            n = (chainlen-i)
            # print n
 
            # calc n
            if n >= len(R2_array):
                R2_array.append(r2)
                n_array.append(2*n-1)
                k_array.append(1)
            else:
                R2_array[n] += r2
                k_array[n] += 1
                n_array[n] = n*2-1
        # print "R = " , R2_array
        # print "n = " , n_array

    R2_array = np.array(R2_array)
    n_array = np.array(n_array)
    k_array = np.array(k_array)
    R2_array /= k_array*n_array
    R2_array = R2_array[:-Noff]
    n_array = n_array[:-Noff]
    save_plot_r2n(n_array, R2_array,psffile,frame,logplot=False)

    return None


    
def calc_r2n(args):
    """
    input : args -> traj from read_traj_vmd
    produces an animation of trajectory relaxation parameter R2(n)/n
    """
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    # Nres = len(u.residues)
    os.system("rm R2*.png")
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        print " frame = %d " % ts.frame
        get_r2n(u,psffile,ts.frame,Noff=args.Noff,logplot=args.logplot)
        # R2 = get_r2n_mono(u,psffile,ts.frame)
    os.system("convert -delay 60 R2%s_*.png animR2%s.gif" % (psffile,psffile))
    os.system("rm R2*.png")

def main():
    
    args = read_traj_vmd()
    calc_r2n(args)
    # print R2


if __name__ == '__main__':
    main()
