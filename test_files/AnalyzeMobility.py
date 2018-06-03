#!/usr/bin/env python

# Program: AnalyzeMobility.py
# Purpose: analyzes rmsd, mobility coeffecient
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeMobility.py for help, 
# Requires: read_parameters.py

# Theory :
# g2(t) mean-square displacement relative to the chains center of mass 
# g2(t) = (1/MN)*(r(t)-r(0) - rc(t) + rc(0))
# D = g2(t) / 6*t when t-> inf 

# The plots seem ok, but they are not what I am expecting 


# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 


import numpy as np
from MDAnalysis import *
from read_parameters import read_traj_vmd
import os
from numba import jit
from save_plots import save_plot_mobility as save_plot





# @jit('nopython=True')
def calc_rmsd(a,b):
    """calculation of difference using numba"""
    N = a.shape[0]
    K = a.shape[1]
    for i in range(N):
      diff = 0.0
      for k in range(K):
        tmp = a[i,k] - b[i,k]
        diff += tmp*tmp
    diff /= float(N)*float(K)
    diff = np.sqrt(diff)
    return diff




def main():
    args = read_traj_vmd()
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    ref_atoms = u.selectAtoms("all")
    traj_atoms = u.selectAtoms("all")
    natoms = traj_atoms.numberOfAtoms()
    frame = []
    g2 = []

  # if performing a mass-weighted alignment/rmsd calculation
  #masses = ref_atoms.masses()
  #weight = masses/numpy.mean(masses)



  # reference centre of mass system
    ref_com = ref_atoms.centerOfMass()
    ref_coordinates = ref_atoms.coordinates() - ref_com
  # diff_coordinates = ref_atoms.coordinates().copy

  # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.coordinates().copy()

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
      # shift coordinates for rotation fitting
      # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass()
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        # diff = calc_rmsd(traj_coordinates,ref_coordinates)
        diff = traj_coordinates - ref_coordinates
        norm = np.linalg.norm(diff,axis=1)
        norm *= norm
      # rms = sqrt(mean_squared_error(traj_coordinates,  ref_coordinates ))
      # diff_coordinates = traj_coordinates - ref_coordinates
      # mean = np.mean(LA.norm(diff_coordinates,axis=0))
      # difference = 
      # print diff_coordinates
      # print "%5d  %8.3f A" % (k, rmsd[k])
        print "frame " , ts.frame, " diff " , np.average(norm)
      # print "frame " , ts.frame, " n " , rms
        frame.append(ts.frame)
        g2.append(np.average(norm))
    frame = np.array(frame)
    g2 = np.array(g2)
    # os.system("convert -delay ")
    save_plot(frame,g2,psffile)  
if __name__ == '__main__':
  main()