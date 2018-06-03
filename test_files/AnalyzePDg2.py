#!/usr/bin/env python

# Program: AnalyzeIndividChain.py
# Purpose: calculates crystallinity of each individual chainlength,
# allowing to analyze polydisperse systems
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeIndividChain.py for help,
# example: python AnalyzeIndividChain.py
# Requires: read_parameters.py

# Theory:
# we have a function that_calculates_everything
# SimData # plots everything,
#   by looping through trajectory and invoking that_calculates_everything
# Simpoints # initializes data for the plot to be imported to the animation
# main:
#   gets psffile,
#   creates a universe
#   creates animation
# if __main__
#   initializes parameters to plot, xlabel and stuff
#   calls main

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import MDAnalysis
import read_parameters
import get_path_names
import AnalyzeChain
import AnalyzeLog


# this part that calculates everything
def analyze_chain_hist(u, args):
    """
    loops through residues, finds crystallinity of each chain, and writes
    an array of chainlengths and chaing2s(chain crystallinities)
    input: Universe (MDAnalysis)
    output: chainlengths, chaing2s
    part of that_calculates_everything
    """
    chainlengths, chaing2s = [], []

    for res in u.residues:
        chords = AnalyzeChain.get_bondlist_coords(res)
        resid_cryst = AnalyzeChain.get_chain_crystallinity(chords)
        # all_g2 += resid_cryst
        chainlengths.append(len(res))
        chaing2s.append(resid_cryst)
    chainlengths = np.array(chainlengths)
    chaing2s = np.array(chaing2s)
    return chainlengths, chaing2s


def get_diff_chain_cryst(u, args):
    """
    creates array chaintypes which is the unique chainlengths.
    Then histogramms that crystallinities
    input: u, args
    output: chaintypes, histogram
    histgoraming is being done this way

    # i - loop through chaintypes, total Ntypes
    # iall - loop through chainlengths, total Nchains
        calculate crystallinity of each class of chaintypes - i
        save to histogram
    part of that_calculates_everything

    """
    # now we need to histogram the data
    chainlengths, chaing2s = analyze_chain_hist(u, args)
    chaintypes, chaincounts = np.unique(chainlengths, return_counts=True)
    # will be used as bins
    histogram = 0.0*chaintypes
    Ntypes = chaintypes.shape[0]
    Nchains = chainlengths.shape[0]
    # iall - loop through chainlengths, total Nchains
    # i - loop through chaintypes, total Ntypes
    for i in range(Ntypes):
        cryst = 0.
        kk = 0
        for iall in range(Nchains):
            if chainlengths[iall] == chaintypes[i]:
                cryst += chaing2s[iall]
                kk += 1
        if kk == 0:
            cryst = 0.
        else:
            cryst /= float(kk)
        histogram[i] = cryst
    histogram = np.asarray(histogram)
    return chaintypes, histogram, chaincounts


# this part that plots everything
def SimData(u, args):
    """
    input: u, args
    plot through trajectory and plot frames
    get chaintype, histogram, chaincounts for current frame
    then we plot a histogram

    """
    name = get_path_names.get_filename(args.psffile)
    total_frames, total_crystallinity = [], []
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        print " frame = %d " % ts.frame
        # chaintype, histogram, chaincounts
        x, y, counts = get_diff_chain_cryst(u, args)
        total_frames.append(ts.frame)
        total_crystallinity.append(y)  # y = histogram
        plt.axhline(y=y.min(),
                    xmin=x.min(),
                    xmax=x.max(),
                    linewidth=0.8,
                    color='k',
                    ls='--')
        plt.bar(x, y,
                width=0.2*(x.max() - x.min())/float(len(x)),
                facecolor='green',
                alpha=0.7,
                align='center')
        yield x, y

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
    total_frames = np.asarray(total_frames, dtype='float32')
    total_crystallinity = np.asarray(total_crystallinity, dtype='float32')
    factor = 1e6
    total_frames *= dumpskip*timestep/float(factor)
    for j in range(total_crystallinity.shape[1]):
        np.savez('indchainc'+str(x[j])+'n'
                 + str(counts[j]) +
                 'PDI1.0'+'extra'+name,
                 total_frames, total_crystallinity[:, j])
    # total_crystallinity = np.asarray(total_crystallinity)


def SimPoints(SimData):
    """
    initializes data for the plot to be imported to the animation
    """
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
                                  blit=True, interval=150,
                                  repeat=True)

    ani.save('Cryst{filename}.mp4'.format(filename=psffile))
    plt.close()

if __name__ == '__main__':
    # here are the global constants (therefore all capital letters)
    FIG = plt.figure()
    AX = FIG.add_subplot(111)
    LINE, = AX.plot([], [], '--')
    AX.set_ylim(0, 1)
    AX.set_xlabel(r'$\mathrm{Chain\ length}$')
    AX.set_ylabel(r'$\mathrm{Crystallinity}$')
    args = read_parameters.read_traj_vmd()
    main(args)
