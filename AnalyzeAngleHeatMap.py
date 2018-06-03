#!/usr/bin/env python

import os
import numpy as np
# from MDAnalysis import *
import MDAnalysis
# from read_parameters import read_traj_vmd
import read_parameters
import save_plots
import termcolor
import AnalyzeChain
import get_path_names
import AnalyzeLog
# import numba
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# this alghoirthm looks for the lamellae detection:
def figsize(scale):
    # Get this from LaTeX using \the\textwidth
    # will make the plot and fonts the right size
    ########################

    # fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
    fig_width_pt = 500.0  # Get this from LaTeX using \showthe\columnwidth
    #########################
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    # golden_mean = 0.8

    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean       # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size



pgf_with_latex = {
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 16,
    "text.fontsize": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "figure.figsize": figsize(1.0),
    "figure.autolayout": True}

# Optionally set font to Computer Modern to avoid common missing font errors
# matplotlib.rc('font', family='serif', serif='cm10')

# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

matplotlib.use('pgf')
matplotlib.rcParams.update(pgf_with_latex)


def main(args):

    #
# plt.subplot(2, 1, 1)
# im = plt.pcolormesh(x, y, z, cmap=cmap, norm=norm)

    plt.clf()
    timeseries = []
    psffile = get_path_names.get_filename(args.psffile)
    u = MDAnalysis.Universe(psffile+'.psf', args.traj)
    res = u.selectAtoms("resid 50")
    nframes = 0
    Nbonds = len(u.residues[0])
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        nframes += 1
        print "frame = {frame}, percent of lamellae={num_of_lam}"\
               .format(frame=ts.frame,
                       num_of_lam=1.)
        # angles_res = np.log(res.angles.angles() + 1e-10)
        angles_res = (180.*res.angles.angles()/np.pi + 1e-10)
        timeseries.append(angles_res)
    timeseries = np.array(timeseries)

    cmap = plt.get_cmap('brg')
    levels = MaxNLocator(nbins=15).tick_values(timeseries.min(), timeseries.max())
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    plt.imshow(timeseries, interpolation='sinc',cmap=cmap, aspect=Nbonds/float(nframes))
    # plt.imshow(timeseries, origin='lower', interpolation='nearest', aspect=Nbonds/float(nframes))
    plt.xlabel('$Monomer\ index$')
    plt.ylabel('$Time, steps$')
    plt.title('Angle distribution')
    plt.colorbar()
    # plt.show()
    plt.savefig("angle"+psffile+".pdf")

if __name__ == '__main__':

    args = read_parameters.read_traj_vmd()
    main(args)
