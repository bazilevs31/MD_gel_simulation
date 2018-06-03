#!/usr/bin/env python

# Program: CalcPolydispIndex.py
# Purpose: plots the histogram with polydispersity index mentioned, from the  *.npz files using Python
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CalcPolydispIndex.py -h for help,
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

from read_parameters import read_calc_pdi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def main():
    """
    uses read_parameters to get: args.npzfile, args.mu, args.sigma
    calculates the PDI parameter and replots the histogram of the PD melt
    with the PDI parameter
    """
    args = read_calc_pdi()
    data = np.load(args.npzfile)
    n = data['arr_1']
    bins = data['arr_0']
    mu,sigma = args.mu, args.sigma
    # n, bins, patches = plt.hist(x, bins=self.Nbins,facecolor='green', alpha=0.75)
    # add a 'best fit' line
    Natoms = (n*bins).sum()
    Nchains = n.sum()
    N1 = 0; N2 = 0; N3 = 0
    for i in range(len(n)):
        N1 += n[i]*bins[i]**2
        N2 += n[i]*bins[i]
        N3 += n[i]

    PDI = (float(N1)/float(N2))/(float(N2)/float(N3))

    print ("current directory is")
    print(os.getcwd() + "\n")
    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    curdir = curdir
    figuresdir = figuresdir

    # fig = plt.figure()

    # width = .35
    print "PDI", PDI

    plt.bar(bins, n, width = 30,facecolor='green', alpha=0.75)

    plt.xlabel(r'$\mathrm{Chain\ length}$')
    plt.ylabel(r'$\mathrm{Number\ of\ chains}$')
    # plt.title(r'$\mathrm{Polydisperse\ chains\ with }\ \mu=%d,\ \sigma=%d,\ PDI=%5.3f,\ N_{atoms}=%d,,\ N_{chains}=%d$' % (mu,sigma,PDI,Natoms,Nchains))
    plt.axis([bins.min()-10, bins.max()+50, 0, n.max() + int(0.05*n.max())])
    plt.grid(True)
    plt.savefig(figuresdir + '/pdchainPDI%5.3f.pdf' % (PDI))
    # np.savez(figuresdir + '/pdmu%ds%dPDI%5.3f' % (mu,sigma,PDI), bins, n)
    # plt.show()

if __name__ == '__main__':
    main()
