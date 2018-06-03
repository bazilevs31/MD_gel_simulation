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
from MDAnalysis import *
from read_parameters import read_traj_vmd
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from AnalyzeChain import get_chain_crystallinity,get_bondlist_coords
from save_plots import save_plot_cryst

from itertools import cycle,izip

# I need to make an animation
# where different chain lengths are being compared for their crystallinity

def plot_hist(bins, histogram, frame, psffile):
    """
    plots histgoram for ndisp melt, can handle even a monodisperse case
    bins = chaintypes,
    histogram = crystallinity obtained from chaing2s.
    this will plot the histogram
    and then convert it to gif
    """
    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    N = len(histogram)
    xmin = 0.8*bins.min()
    xmax = 1.2*bins.max()
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    plt.ylim(0,1)
    plt.xlim(xmin,xmax)
    plt.xlabel(r'$\mathrm{Chain\ length}$')
    plt.ylabel(r'$\mathrm{crystallinity}$')
    plt.grid(True)
    plt.title('%s chain crystallinity time=%f ' % (psffile,frame))
    plt.axhline(y=histogram.min(), xmin=xmin-30, xmax=xmax, linewidth=0.8, color = 'k',ls='--')
    plt.bar(bins,histogram, width=0.2*(xmax-xmin)/float(N), facecolor='green', alpha=0.7,align='center')
    plt.legend(loc='best')
    plt.savefig('hist_%.5d.png' % int(frame))
    plt.cla()


def get_hist_cryst(chainlengths,chaing2s, frame,psffile):
    """
    create  the same thing with unique

    chaintypes : unique elements of chainlengths
    chaing2s: chain crystallinities
    frame : current time frame
    psffile : file to write animation to
    # i - loop through chaintypes, total Ntypes
        # iall - loop through chainlengths, total Nchains
            calculate crystallinity of each class of chaintypes - i
            save to histogram
    use plot_hist to plot the results
    """
    chaintypes = np.unique(chainlengths)  # will be used as bins
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
        if kk==0:
            cryst = 0.
        else:
            cryst /= float(kk)
        histogram[i] = cryst
    plot_hist(chaintypes,histogram,frame,psffile)
    return chaintypes, histogram


def get_diff_chain_cryst(chainlengths,chaing2s, frame,psffile):
    """
    this function analyzes crystallinities of all chains
    creates an intermediate array
    C_array of size (Nframes, N chaintypes + 1)
    frame chaintype=1 chaintype=2 ....
     |      |              |

    histgoraming is being done this way

    # i - loop through chaintypes, total Ntypes
    # iall - loop through chainlengths, total Nchains
        calculate crystallinity of each class of chaintypes - i
        save to histogram

    """
    chaintypes = np.unique(chainlengths) # will be used as bins
    histogram = 0.0*chaintypes
    frames = []
    Ntypes = chaintypes.shape[0]
    Nchains = chainlengths.shape[0]
    # iall - loop through chainlengths, total Nchains
    # i - loop through chaintypes, total Ntypes
    frames.append(frame)
    for i in range(Ntypes):
        cryst = 0.
        kk = 0
        for iall in range(Nchains):
            if chainlengths[iall]==chaintypes[i]:
                cryst += chaing2s[iall]
                kk += 1
        if kk==0:
            cryst = 0.
        else:
            cryst /= float(kk)
        histogram[i] = cryst
        C_array.append(histogram)
    C_array = np.array(C_array)
    frames = np.array(frames)
    plt.xlim(0.,1.)
    for j in range(chaintypes):
        plt.plot(frames[:],C_array[:],label=str(chaintypes[j]))
    return None

def analyze_chain_hist(args):
    """
    function that does everything
    """
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    Nres = len(u.residues)
    print "number of residues is %r" % Nres
    frames = []
    C_array = []
    g2 = []
    chainlengths = []
    chaing2s = []
    res_g2 = 0.



    os.system("rm hist_*")
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        all_g2 = 0.
        k = 0
        chainlengths = []
        chaing2s = []


        frames.append(ts.frame)
        for res in u.residues:
            chords = get_bondlist_coords(res)
            res_g2 = get_chain_crystallinity(chords, threshold=args.threshold, neigh=args.neigh)
            all_g2 += res_g2
            chainlengths.append(len(res))
            chaing2s.append(res_g2)
        chainlengths = np.array(chainlengths)
        chaing2s = np.array(chaing2s)
        bins, histogram = get_hist_cryst(chainlengths,chaing2s,ts.frame,psffile)
        C_array.append(histogram)

        print ("frame %d , g2 = %f" %(ts.frame, all_g2 / float(Nres)))
        g2.append(all_g2 / float(Nres))

    C_array = np.array(C_array)
    frames = np.array(frames)
    plt.ylim(0.,1.)
    plt.xlabel(r'$\mathrm{time}$')
    plt.ylabel(r'$\mathrm{crystallinity}$')
    plt.grid(True)
    plt.title(' chain crystallinity evolution ')
    lines = [u'D', u'o', u'^', u'>', u's', u'8', u'<', u'>',  u'*',  u'H', u'h', u'p', u'v', u'D', u'd',"-","--","-.",":"]
    colors=['k','y','m','c','b','g','r','#aaaaaa']
    # colors=['m','c','b','g','m','k']

    linecycler = cycle(lines)
    colorcycler = cycle(colors)

    for j in range(len(bins)):
        # plt.plot(frames[:],frames[:])
        clr=next(colorcycler)

        name = 'N = % r' % (bins[j])
        plt.plot(frames[:], C_array[:,j],next(linecycler),color=clr,linewidth=1.2,label=name)
        np.savez('ChainSelection' + psffile + name, frames[:],C_array[:,j])

        # print bins[j]
        # plt.plot(frames[:], C_array[:],label=str(j))
    plt.legend(loc='best')
    plt.savefig('allchain_diff.pdf')
    np.savez('lasthist' + psffile, bins, histogram)
    frames = np.array(frames)
    g2 = np.array(g2)
    os.system("convert -delay 50 hist*.png hist%s.gif" % psffile)
    return None


def main():
    args = read_traj_vmd()
    analyze_chain_hist(args)

if __name__ == '__main__':
    main()
