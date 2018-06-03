#!/usr/bin/env python

# Program: AnalyzeSq.py
# Purpose: analyzes static structure factor
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeSq.py -h for help,
# example: python AnalyzeSq.py -f melt.psf -t melt.dcd
# Requires: read_parameters.py, save_plots

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


# how to test?
# make a lj fluid, -> a lot of atoms
# get their g(r) using VMD
# calculate s(q) using g(r)
# then calculate it using this -> then compare

import numpy as np
from MDAnalysis import *
from read_parameters import read_traj_vmd
import os
import save_plots
from numpy.fft import fftn, fftshift
import scipy
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import get_path_names
# from AnalyzeLog import get_dumpskip, get_timestep/



def calc_sq(args):
    """
    input: args obtained from read_traj_vmd
    calculates sq static parameter to analyze trajectories
    output:
    """
    # psffile = os.path.splitext(args.psffile)[0]
    u = Universe(args.psffile, args.traj)
    psffile = get_path_names.get_filename(u.filename)

    # Ndiv = args.nsub
    Natoms = len(u.atoms)
    # Ndiv = int(Natoms**1/3.)
    Ndiv = 201
# create arrays of positions, get the box dimenstions
    u.trajectory[args.startframe]
    u.SYSTEM.packIntoBox()
    X = u.atoms.positions
    box = u.trajectory.ts.dimensions[:-3]
    length_x = box[-1]  # be careful here
# x = np.arange(-length_x,length_x,dx)
    x = np.linspace(-length_x, length_x, Ndiv+1, endpoint=True)

    f, edges = np.histogramdd(X, bins=(Ndiv, Ndiv, Ndiv))
    delta = x[1]-x[0]

# fast fourier transform, go to reciporal space

    ftk = (fftshift(fftn(fftshift(f))*delta))
    # sk = np.abs(ftk**2) / float(Natoms)
    sk = np.abs(ftk**2) / float(Natoms)
# / float(ndata) probably number of atoms

    # plt.cla()

# basis in reciporal space
    omega = 2*np.pi*np.arange(Ndiv-1) / (length_x)
    omega -= omega[int(Ndiv/2)-1]
    k1, k2, k3 = np.meshgrid(omega, omega, omega)
    kmax = np.sqrt(3.)*omega[-1]
    # print "dk =" , dk, " kmax = ", kmax
    C = norm_sq(sk, k1, k2, k3, Ndiv, kmax)
# we are interested in spherically symmetric case,
# i.e only the module of [k1,k2,k3] = normk is important
    # print C[:,1]
    # plt.plot(C[:,1], C[:,2])
    Nfirst = 10
    save_plots.save_plot(C[Nfirst:-Nfirst, 1],
                         C[Nfirst:-Nfirst, 2],
                         xlabel='q, 1/sigma',
                         ylabel='sq',
                         title='static structure factor',
                         plotname='sq',
                         name=psffile+'f'+str(args.startframe),
                         ymax=1.0, ymin=0.)
    return None


def norm_sq(sk, k1, k2, k3, Ndiv, kmax):
    """
    calculates normk = norm of each vector in grid=|k1,k2,k3|
    bins in, calculating average sq of each bin
    output: C[Nbins,3] = binindex, kval = normk, sq
    """
    # Nbins = int(Ndiv*1)
    normk = np.sqrt(k1*k1 + k2*k2 + k3*k3)
    # array to histogram - dk,kmax,Nbins
    Nbins = int(1.5*Ndiv)
    dk = kmax/float(Nbins)
    C = np.ones((Nbins+1, 3))
    for i in range(Ndiv-1):
        for j in range(Ndiv-1):
            for k in range(Ndiv-1):
                kval = normk[i, j, k]
                bindex = int(kval/dk)
                C[bindex, 0] += 1
                C[bindex, 1] += kval
                C[bindex, 2] += sk[i, j, k]
    C[:, 2] /= C[:, 0]
    C[:, 1] /= C[:, 0]
    # plt.plot(C[:,1], C[:,2])
    # plt.show()
    # print C[:,1]
    # print C[:,2]
    # plt.savefig('lala.pdf')
    return C


def main():
    args = read_traj_vmd()
    # calc_sq_rdf(args)
    # print "doing sq direct"
    calc_sq(args)

if __name__ == '__main__':
    main()
