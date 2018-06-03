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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from read_parameters import read_plo
import read_parameters
import get_path_names
# def del_end(x):
#     """
#     delete npz dimension of the file
#     """
#     return os.path.splitext(x)[0]
# def del_path(x):
# 	return os.path.basename(x)

# def is_valid_file(parser, arg):
#     """
#     Check if arg is a valid file
#     """
#     argpath = os.path.abspath(arg)
#     if not os.path.exists(argpath):
#         parser.error("The file %s doesn't exist " % argpath)
#     else:
#         return arg
# if args.bars:
#         plot_bars(args)
#     elif args.points:
#         tcryst_array = []

#                 flag_cryst = False
#         for jj in range(len(p)):
#             if flag_cryst==False:
#                 if p[jj]>0.3:
#                     flag_cryst=True
#                     tcryst_array.append(t[jj])
#             if args.boxes:
#             bbox_props = dict(boxstyle="round", fc=clr, ec="0.01", alpha=0.2)
#             plt.text(t[-10]-10.0, p[-10]+0.03, j, ha="center", va="center", rotation=15, size=10, bbox=bbox_props)
#             plt.text(t[place], p[place]+yoffset, j, ha="center", va="center", rotation=15, size=10, bbox=bbox_props)
#         if place>40:
#             yoffset = 0.03
#             place = 20
#     # plt.plot((0, t[-1]), (0.2, 0.2), 'k--', alpha=0.6)


#         # plt.cla()

#         # x = filenames
#         # print "x = ", x
#         # print "t = ", tcryst_array
#         # tcryst_array = np.array(tcryst_array)
#         # plt.xlabel('$\mathrm{M_w}$')
#         # plt.ylabel('${0[xlabel]}$'.format(parameters))
#         # plt.title('${0[title]}$'.format(parameters))
#         # plt.plot(tcryst_array,'go-')
#         # plt.xticks(range(len(x)), x)
#         # plt.grid()
#         # # plt.show()
#         # plt.savefig('indtime{0[title]}.pdf'.format(parameters))


#         # plt.bar(filenames,tcryst_array)


def plot_bars(args):
    """plots different histgrams"""
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

    plt.ylabel('${0[ylabel]}$'.format(parameters))
    plt.xlabel('${0[xlabel]}$'.format(parameters))
    plt.title('${0[title]}$'.format(parameters))

    widths = [4, 6, 8, 10, 12]
    colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']

    colorcycler = itertools.cycle(colors)
    widthscycler = itertools.cycle(widths)

    for (i, j) in itertools.izip(args.collection, filenames):
        npz = np.load(i)
        t = npz[args.name_x]
        p = npz[args.name_y]
        clr = next(colorcycler)
        w = next(widthscycler)
        plt.bar(t, p, width=w, facecolor=clr,
                alpha=0.4, label=j, align='center')
    if args.legend:
        plt.legend(loc='best')
    plt.savefig(args.titlename+".pdf")
    return None


def plot_points(args):
    """
    plots stuff with points
    #algorithm for induction time.
    filenames will represent the mol weight, and will be the X axis
    the flag = False will represent that the crystallization has started
    the tcryst = will be the time of crystallization, it will be appended
    to tcryst_array list
    """
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

    lines = [u'D-', u'o-', u'^-', u'>-', u's', u'8', u'<', u'>',  u'*',
             u'H', u'h', u'p', u'v', u'D', u'd', "-", "--", "-.", ":"]
    # widths = [4, 6, 8, 10, 12]
    colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']

    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)

    plt.figure()
    plt.ylabel('${0[ylabel]}$'.format(parameters))
    plt.xlabel(r'${0[xlabel]}$'.format(parameters))
    # plt.title('${0[title]}$'.format(parameters))
    # plt.ylim(0,1)
    # plt.xlim(10,210)
    for (i, j) in itertools.izip(args.collection, filenames):
        # ChainLength = int(re.search('(\d+)n(\d+)', j).group(1))
        # NumOfChains = int(re.search('(\d+)n(\d+)', j).group(2))
        print j
        # Extra_string = (re.findall('fit(\S+)', j))
        npz = np.load(i)
        t = npz[args.name_x]
        p = npz[args.name_y]
        # p_err = npz[args.name_z]
        p_err = t*0.
        # p = p - p[0]
        clr = next(colorcycler)
        # namefile=get_path_names.get_filename(j)
        # namefile = Extra_string[0]
        namefile = j
        plt.errorbar(t, p,
                     yerr=p_err,
                     color=clr ,
                     fmt='o-',
                     linewidth=1.4,
                     label=namefile)
        # plt.plot(t, p, 'go-',color=clr,linewidth=1.4,label=namefile)
        # plt.plot(t, p, next(linecycler),
        #          color=clr,
        #          linewidth=1.2,
        #          label=namefile)
                 #,
                 # label=r'$M_n = %d,\ N_{chains} = %d$'
                 # % (ChainLength, NumOfChains))
    # if args.legend:
    plt.legend(loc='best')
    plt.savefig('{0[title]}.pdf'.format(parameters))
    return None


def main():
    """
    get files, filenames(without .npz), name_x,name_y, titlename
    list a long file of lines(which represent different styles)
    the same with colors
    then use the itertools and plot them, switch to different style by next()
    """

    args = read_parameters.read_plot()
    print args

    if args.bars:
        plot_bars(args)
    elif args.points:
        plot_points(args)
    else:
        raise ValueError('no input name specified')


if __name__ == '__main__':
    main()
