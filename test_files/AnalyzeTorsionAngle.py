#!/usr/bin/env python

# Program: AnalyzeTorsionAngle.py
# Purpose: analyzes population of each torision angle
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeTorsionAngle.py for help,
# example: python AnalyzeTorsionAngle.py -mu 100 -s 20 -nch 100
# Requires: read_parameters.py, save_plots

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
from MDAnalysis import *
from read_parameters import read_traj_vmd
import os
from save_plots import save_plot_cryst as save_plot
from AnalyzeLog import get_dumpskip, get_timestep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    """main programm"""
    args = read_traj_vmd()
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    # u.trajectory[-1]
    angles = u.angles.angles()
    ymax = len(angles)/10
    xmax = angles.max()*1.1
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        print " frame = %d " % ts.frame
        angles = u.angles.angles()
        plt.cla()
        plt.ylim(0.,ymax)
        plt.xlim(0.,xmax)
        plt.title("frame = %d" % ts.frame)
        plt.hist(angles,bins=100,alpha=0.75)
        plt.savefig("angles_hist_%.5d.png" % ts.frame)

    os.system("convert -delay 60 angles_*.png anim_angles.gif")
    os.system("rm angles*.png")


if __name__ == '__main__':
    main()
