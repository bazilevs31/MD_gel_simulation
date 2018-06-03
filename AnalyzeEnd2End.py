#!/usr/bin/env python

# Program: AnalyzeEnd2End.py
# Purpose: Calculates end-to-end distance, can handle polydisperse systems
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeEnd2End.py -h for help,
# Requires: read_parameters.py

# theory max_chain, get_r2n - are the working horses
# after its all vizualization
# todo : check whether get_r2n works properly

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import MDAnalysis
import read_parameters
import get_path_names


def max_chain(u):
    """
    input: MDAnalysis universe
    output: maximum length of all of the chains present(integer)
    """
    maxlen = 0
    for res in u.residues:
        reslen = len(res)
        if maxlen < reslen:
            maxlen = reslen
    # print "maxlen = %f" % maxlen
    return maxlen - 1


def get_r2n(u):
    """
    create a list of
    R2_array - array of distances
    n_array - array of number of bonds between atoms
    k_array - array of number of atoms with this bonds
    start looping in residues
    for every residue:
        start from the middle of it
            calculate the closest atom_i - atom_-i
            if it is the first time we have the number of bonds so big,
             we expand our lists by appending
            else: we just put it to the nth position
    then the last elements of the array will be deleted,
     since there is not enough statistics for this chains
    """
    Noff = 1
    R2_array = []
    n_array = []
    k_array = []

    for res in u.residues:
        chainlen = len(res)/2
        for i in range(chainlen)[::-1]:
            ag1, ag2 = res.atoms[i].pos, res.atoms[-i-1].pos
            tmpdiff = ag1-ag2
            r2 = np.dot(tmpdiff, tmpdiff)
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
    # save_plot_r2n(n_array, R2_array,psffile,frame,logplot=False)
    return n_array, R2_array


def SimData(u, args):
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        print " frame = %d " % ts.frame
        x, y = get_r2n(u)
        plt.plot(x, y)
        yield x, y


def SimPoints(SimData):
    x, y = SimData[0], SimData[1]
    LINE.set_data(x, y)
    return LINE,


def main(args):
    """
    main func
    """
    u = MDAnalysis.Universe(args.psffile, args.traj)
    psffile = get_path_names.get_filename(args.psffile)

    # Now call the animation package:
    # (simData is the user function
    # serving as the argument for simPoints):
    ani = animation.FuncAnimation(FIG, SimPoints, SimData(u, args),
                                  blit=True, interval=80,
                                  repeat=True)

    ani.save('{filename}.mp4'.format(filename=psffile))
    plt.close()

if __name__ == '__main__':
    # here are the global constants (therefore all capital letters)
    FIG = plt.figure()
    AX = FIG.add_subplot(111)
    LINE, = AX.plot([], [], 'bo')
    AX.set_ylim(0, 40)
    args = read_parameters.read_traj_vmd()
    main(args)
