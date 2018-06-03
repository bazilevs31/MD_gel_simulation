#!/usr/bin/env python

# Program: AnalyzeRgRe.py
# Purpose: calculates Rg, Re parameters
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeRgRe.py -h for help, 
# Requires: read_parameters.py, save_plots

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 

import numpy as np
from MDAnalysis import *
import MDAnalysis.analysis.align as MDrmsd
from read_parameters import read_traj_vmd
import os
from save_plots import save_plot_re, save_plot_rg,save_plot_rmsd


def distance_sq_coords(coord1, coord2):
    diff = np.array(coord1)-np.array(coord2)
    x = np.dot(diff, diff)
    return x

def calc_rmsd(args):
    """
    input : args -> from read_traj_vmd
    calculates g2(t) rmsd, and plots to file using save_plot_rmsd
    """

    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    ref_atoms = u.selectAtoms("all")
    traj_atoms = u.selectAtoms("all")
    frame = []
    rmsd_array = []


  # reference centre of mass system
    ref_com = ref_atoms.centerOfMass()
    ref_coordinates = ref_atoms.coordinates() - ref_com
  # diff_coordinates = ref_atoms.coordinates().copy

  # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.coordinates().copy()

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
      # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass()
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        rmsd_value = MDrmsd.rmsd(traj_coordinates,ref_coordinates)

        print "frame " , ts.frame, " rmsd " , rmsd_value
        frame.append(ts.frame)
        rmsd_array.append(rmsd_value)
    frame = np.array(frame)
    rmsd_array = np.array(rmsd_array)
    save_plot_rmsd(frame,rmsd_array,psffile)
    return None

def calc_rgre(args):
    """
    input : args -> from read_traj_vmd
    calculates Rg(t), Ree(t) and plots it to file using save_plot_re,save_plot_rg
    """
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    reslen = len(u.residues[0])
    frame = []
    Rg_array = []
    Re_array = []

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        Re = np.sum(np.array([distance_sq_coords(myres.atoms[0].pos,myres.atoms[-1].pos) for myres in u.residues]))
        Re /= len(u.residues)
        Rg = np.sum(np.array([myres.atoms.radiusOfGyration(pbc=False) for myres in u.residues]))
        Rg /= len(u.residues)
        frame.append(ts.frame)
        print "(R(last)-R(first))**2 = %f" % Re
        Re /= (reslen-1.) 
        # print ("frame %d , Re^2/(reslen.*Rg^2) = %f" %(ts.frame, Re/(Rg**2)))
        Rg_array.append(Rg)
        Re_array.append(Re)
        print ("frame %d , Rg = %f, R2e/N = %f, R2ee/Rg**2 =  %f, c_inf = (Re^2/n)/(l**2) = %f" %(ts.frame, Rg, Re, (reslen-1.)*Re/Rg**2, Re/(0.85**2)))

    frame = np.array(frame)
    Rg_array = np.array(Rg_array)
    Re_array = np.array(Re_array)
    save_plot_re(frame, Re_array,psffile)
    save_plot_rg(frame, Rg_array,psffile)
    return None

def main():

    args = read_traj_vmd()
    calc_rgre(args)
    calc_rmsd(args)



if __name__ == '__main__':
    main()