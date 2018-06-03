#!/usr/bin/env python

# Program: AnalyzeSq.py
# Purpose: analyzes static structure factor
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeSq.py -h for help,
# example: python AnalyzeSq.py -
# Requires: read_parameters.py, save_plots

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
from MDAnalysis import *
from read_parameters import read_traj_vmd
import os
from save_plots import save_plot_vs_time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy.fft import fftn, fftshift
# from AnalyzeLog import get_dumpskip, get_timestep/

# args = read_traj_vmd()
# print args
def calc_sq(args):
    """
    input: args obtained from read_traj_vmd
    calculates sq static parameter to analyze trajectories
    output:
    """
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    Ndiv = args.nsub
    # psffile = "poly_40.psf"
    # traj = "wraptraj_40.dcd"
    u.trajectory[-1]
    # Nsub = args.nsub
    # cryst_all = 0.0
    X = u.atoms.positions
    box = u.trajectory.ts.dimensions[:-3]
    length_x = box[-1] # be careful here
    x=np.linspace(-length_x,length_x,Ndiv+1,endpoint=True)
    delta=x[1]-x[0]
    y=np.linspace(-length_x+delta/2.,length_x-delta/2.,Ndiv+1,endpoint=True)
    # H, edges = np.histogramdd(X, bins = (Ndiv, Ndiv, Ndiv))
    f, edges = np.histogramdd(X, bins = (Ndiv, Ndiv, Ndiv))

    L2 = np.abs(y[-1]-y[0])
    omega = 2*np.pi*np.arange(Ndiv-1)/(L2+delta)

    omega -= omega[int(Ndiv/2)-1]

    ftk = (fftshift(fftn(fftshift(f))*delta))
    sk = np.abs(ftk**2) # / float(ndata) probably number of atoms


    # plt.cla()

    k1,k2,k3= np.meshgrid(omega,omega,omega)
    normk = np.sqrt(k1*k1 + k2*k2 + k3*k3)
    print normk
    Nbins = Ndiv*1.5
    kmax = np.sqrt(3)*omega[-1]
    dk = kmax/Nbins
    print Nbins
    # % Exclude points where full sphere of radius normk not within box
    #  Radial Average
       # % Bin in k-space
    C=np.ones((Nbins+2,3))

    for i in range(Ndiv-1):
        for j in range(Ndiv-1):
            for k in range(Ndiv-1):
                kval=normk[i,j,k]
                bindex = int(kval/dk)
                if bindex>Nbins+2:
                    continue
                C[bindex,0] += 1
                C[bindex,1] += kval
                C[bindex,2] += sk[i,j,k]
    C[:,1] /= C[:,0]

    plt.plot(C[:,1], C[:,2])
    save_plot_vs_time(C[:,1], C[:,2], "S(q) structure factor %s" % psffile, xtitle="q,1/distance",ytitle="S(q)")
    return None


def main():
    args = read_traj_vmd()
    calc_sq(args)

if __name__ == '__main__':
    main()
