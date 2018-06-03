#!/usr/bin/env python

# Program: PlotAll.py
# Purpose: plot various *.npz files using Python
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python PlotAll.py for help,
# example: python PlotAll.py --points -a Plot1.npz -a Plot2.npz -t <title>
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import os
# from itertools import itertools.cycle, itertools.izip
import itertools
import re

import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from read_parameters import read_plo
import read_parameters
import get_path_names

fig_width_pt = 300 # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width, fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 10,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': True,
          'figure.figsize': fig_size,
          "figure.autolayout": True}
plt.rcParams.update(params)


# def figsize(scale):
#     # Get this from LaTeX using \the\textwidth
#     # will make the plot and fonts the right size
#     ########################

#     # fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#     fig_width_pt = 350.0  # Get this from LaTeX using \showthe\columnwidth
#     #########################
#     inches_per_pt = 1.0/72.27               # Convert pt to inches
#     golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
#     # golden_mean = 0.8

#     fig_width = fig_width_pt*inches_per_pt  # width in inches
#     fig_height = fig_width*golden_mean       # height in inches
#     fig_size = [fig_width, fig_height]
#     return fig_size


#    3 fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#    4 inches_per_pt = 1.0/72.27               # Convert pt to inch
#    5 golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
#    6 fig_width = fig_width_pt*inches_per_pt  # width in inches
#    7 fig_height = fig_width*golden_mean      # height in inches
#    8 fig_size =  [fig_width,fig_height]
#    9 params = {'backend': 'ps',
#   10           'axes.labelsize': 10,
#   11           'text.fontsize': 10,
#   12           'legend.fontsize': 10,
#   13           'xtick.labelsize': 8,
#   14           'ytick.labelsize': 8,
#   15           'text.usetex': True,
#   16           'figure.figsize': fig_size}
#   17 pylab.rcParams.update(params)


# pgf_with_latex = {
#     "pgf.texsystem": "pdflatex",
#     "text.usetex": True,
#     "font.family": "serif",
#     "axes.labelsize": 16,
#     "text.fontsize": 12,
#     "legend.fontsize": 12,
#     "xtick.labelsize": 16,
#     "ytick.labelsize": 16,
#     "figure.figsize": figsize(1.0),
#     "figure.autolayout": True}

# Optionally set font to Computer Modern to avoid common missing font errors
# matplotlib.rc('font', family='serif', serif='cm10')

# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

# matplotlib.use('pgf')
# matplotlib.rcParams.update(pgf_with_latex)


args = read_parameters.read_plot()

# def plot_points(args):
#     """
#     plots stuff with points
#     #algorithm for induction time.
#     filenames will represent the mol weight, and will be the X axis
#     the flag = False will represent that the crystallization has started
#     the tcryst = will be the time of crystallization, it will be appended
#     to tcryst_array list
#     """
filenames = [get_path_names.get_filename(x) for x in args.collection]
parameters = {'title': args.titlename,
              'xlabel': args.xlabel,
              'ylabel': args.ylabel,
              'dumpskip': args.dumpskip}

if args.logplot:
    plt.xscale('log')
    parameters['xlabel'] = ''.join(('log', 'parameters[\'xlabel\']'))
    print "xlogscale"
else:
    print "regular xscale"

lines = [u'D-', u'o-', u'^-', u'>-', u's-', u'8-', u'<-', u'>',  u'*',
         u'H', u'h', u'p', u'v', u'D', u'd', "-", "--", "-.", ":"]
# widths = [4, 6, 8, 10, 12]
colors = ['yellow','saddlebrown','green','blue', 'k','c', 'g', 'r','b','#FF8C00', 'y', 'm',  '#aaaaaa']

linecycler = itertools.cycle(lines)
colorcycler = itertools.cycle(colors)

#     fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.semilogx(data[:,1],data[:,2])

# ax2 = plt.axes([.65, .6, .2, .2], axisbg='y')
# ax2.semilogx(data[3:8,1],data[3:8,2])
# plt.setp(ax2, xticks=[], yticks=[])

plt.figure()
# plt.ylabel(r'${0[ylabel]}$'.format(parameters))
plt.ylabel(r'${0[ylabel]}$'.format(parameters))
plt.xlabel(r'${0[xlabel]}$'.format(parameters))
# plt.title(r'${0[title]}$'.format(parameters))
plt.ylim(0.0,0.6)
plt.xlim(0,0.3)
# plt.xlim(0,0.3)
# plt.gca().tight_layout()
# plt.xlim(1)
xloc = 0.00

k = 0

for (i, j) in itertools.izip(args.collection, filenames):



    # ChainLength = int(re.search('(\d+)n(\d+)', j).group(1))
    # NumOfChains = int(re.search('(\d+)n(\d+)', j).group(2))
    print j
    # Extra_string = (re.findall('fit(\S+)n_avrami', j))
    npz = np.load(i)
    t = npz[args.name_x]
    p = npz[args.name_y]
    # p_err = npz[args.name_z]
    p_err = t*0.
    # p = p - p[0]
    mrk = next(linecycler)
    clr = next(colorcycler)
    # xloc += 0.05
    # yloc = np.where(abs(t-xloc)<0.01)
    # if k == 0:
    #   namefile = "monodisperse, $M_w$=40"
    #   padd = p[yloc[0]][0]+0.08
    # elif k==1:
    #   namefile = "bidisperse $(1)$"
    #   padd = p[yloc[0]][0]+0.04
    # elif k==2:
    #   namefile = "bidisperse $(2)$"
    #   padd = p[yloc[0]][0]+0.02
    # elif k==3:
    #   namefile = "monodisperse, $M_w$=160"
    #   padd = p[yloc[0]][0]+0.02
    # else:
    #   print "error"

    # ChainLength = int(re.search('(\d+)n(\d+)', j).group(1))
    # ChainLength = (re.search('fit(\S+)chi0', j).group(1))
    # print ChainLength
    # namefile=get_path_names.get_filename(j)
    # namefile = str(ChainLength)
    namefile = j
    # if namefile in ['x1','x0.5','x10']:
    #   namefile = "cooling $"+namefile+"$"

    # if ChainLength>100:
    # if k>1:
    # padd = -0.04
    # if k>1:
    # if k>1:
    #   xloc += 0.7
    # else:
    #   xloc += 0.11
    # namefile = Extra_string[0]
    # namefile = j
    print t
    print p
    # plt.xlim(0.3, 0.8)
    plt.errorbar(t, p,
                 yerr=p_err,
                 color=clr,
                 fmt=mrk,
                 linewidth=1.8,
                 label=namefile)



    # plt.markevery(5)

    # print yloc
    # print yloc[0]
    # print p[yloc[0]]
    # plt.plot(t, p, mrk,color=clr,markevery=1,linewidth=1.5,label=namefile)
    # plt.annotate(namefile, xy=(0.27, p.max()), xytext=(0.27, 1.2*p.max()))
    # plt.text(xloc, padd, namefile, style='italic',
        # bbox={'facecolor':clr, 'alpha':0.5, 'pad':7})
    # k+=1

    # plt.plot(t, p, next(linecycler),
    #          color=clr,
    #          linewidth=1.2,
    #          label=namefile)
             #,
             # label=r'$M_n = %d,\ N_{chains} = %d$'
             # % (ChainLength, NumOfChains))
# if args.legend:
# plt.legend(loc='best')
# plt.legend(loc='upper center')
plt.savefig('{0[title]}.pdf'.format(parameters))
# return None

