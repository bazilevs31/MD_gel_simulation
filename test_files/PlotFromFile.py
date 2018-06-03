#!/usr/bin/env python

# Program: PlotFromFile.py
# Purpose: plot various *.npz files using Python
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python PlotFromFile.py for help,
# example: python PlotFromFile.py  -in files.input -p files.param
# -bar (True False)
# files.input:
# name of .npz file -> final name in the plot
# files.params:
# xlabel time...
# ylabel crystallinity
# Requires: read_parameters.py


# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


# todo : adding the maximum value of the plot, as additional parameter
# also it would be go if there was only one .param file, and the output file
# (the only difference in different .param files) would be speciafied via -t
# key

# todo
# to handle the parsing through chain info more properly
#

# lala = {
#     "something": "something"
#     "chain_lengths": [1,2,3,...]
#     "numchains" : [13,3,]
# }

import os
import itertools
import re
import numpy as np

import termcolor
import read_parameters
import get_path_names
import scipy.optimize
import ast


import matplotlib
matplotlib.use('Agg')

publication_quality = True
if publication_quality:
    print "publication_quality"
    def figsize(scale):
        # Get this from LaTeX using \the\textwidth
        # will make the plot and fonts the right size
        ########################

        fig_width_pt = 350.  # Get this from LaTeX using \showthe\columnwidth
        #########################
        inches_per_pt = 1.0/72.27               # Convert pt to inches
        golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio

        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*golden_mean       # height in inches
        fig_size = [fig_width, fig_height]
        return fig_size

    pgf_with_latex = {
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,
        "font.family": "serif",
        "axes.labelsize": 12,
        "text.fontsize": 12,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.figsize": figsize(1.0),
        "figure.autolayout": True}

    matplotlib.use('pgf')
    matplotlib.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

def func(x, x0, kappa, n):
    """
    input: x (array), x0,kappa, n - parameters
    output: func(x)
    """
    return x0 * (1. - np.exp(-kappa * x**n))


def read_plot_parameters(args):
    """
    reads a text files with two entries per line
    creates a dictionary from file.
    input: args (a dict obtained from read_parameters module)
    output: parameters (dict)

    parameters: a dict of plot parameters that ready to use

    Example:

    {'normalize': False,
    'figuresdir': '/Users/bazilevs/Desktop/NewSimulation/results_march31/
     npz_crystllinities/figures',
    'ymax': 0.57,
    'fit': True,
    'title': 'newtest',
    'xlog': False,
    'plotname': 'ICC',
    'ylog': False,
    'xlabel': 'time,10^6lj',
    'ylabel': 'chi',
    'ymin': 0.0,
    'legend': True,
    'parse_chain_info': True,
    'dumpskip': 15000,
    'timestep': 0.005,
    'scale_dump': True,
    'plot_style': lines/points/linespoints
    'name': 'newtest'}
    """
    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    filename = args.param
    with open(filename, 'r') as document:
        tmp_parameters_dict = {}
        for line in document:
            if line.strip():  # non-empty line?
                # None means 'all whitespace', the default
                key, value = line.split(None, 1)

                if key in ['normalize', 'fit',
                           'xlog', 'ylog',
                           'parse_chain_info',
                           'scale_dump']:
                    tmp_parameters_dict[key] = ast.literal_eval(value)
                    # tmp_parameters_dict[key] = bool(value)
                else:
                    tmp = value.split()[0]  # first element of the rest
                    tmp_parameters_dict[key] = tmp
                    if key == 'ymin':
                        tmp_parameters_dict[key] = float(tmp)
                    elif key == 'ymax':
                        tmp_parameters_dict[key] = float(tmp)

    parameters = dict()
    parameters['plotname'] = tmp_parameters_dict.get('plotname', 'ICC')
    parameters['xlabel'] = tmp_parameters_dict.get('xlabel',
                                                   'time, 10^6 lj units')
    parameters['ylabel'] = tmp_parameters_dict.get('ylabel', 'crystallinity')
    parameters['xlog'] = tmp_parameters_dict.get('xlog', False)
    parameters['ylog'] = tmp_parameters_dict.get('ylog', False)
    parameters['legend'] = tmp_parameters_dict.get('legend', True)
    parameters['plot_style'] = tmp_parameters_dict.get('plot_style',
                                                       'linespoints')
    parameters['scale_dump'] = tmp_parameters_dict.get('scale_dump', False)
    parameters['timestep'] = tmp_parameters_dict.get('timestep', 0.005)
    parameters['dumpskip'] = tmp_parameters_dict.get('dumpskip', 15000)
    parameters['ymin'] = tmp_parameters_dict.get('ymin', 0.0)
    parameters['ymax'] = tmp_parameters_dict.get('ymax', 1.0)
    parameters['figuresdir'] = figuresdir

    parameters['title'] = args.titlename
    parameters['name'] = args.titlename

    # divide by x_0 or no?
    parameters['normalize'] = tmp_parameters_dict.get('normalize', False)
    parameters['fit'] = tmp_parameters_dict.get('fit', False)
    parameters['parse_chain_info'] = tmp_parameters_dict.get(
        'parse_chain_info', False)
    print termcolor.colored(parameters, 'blue')
    return parameters


def parse_file_info(filename):
    """
    parses through chain length info from file.in file_in

    input: filename: string, parse_chain_info: bool
    output: all_files_info_list: dict

    reads all the files required to plot
    parses through the string and returns a dict which has
    name - bi, pd, md
    all_files_info_list = {'chains_info':[(40,1500),(160,1500)],
    'extrastring': 'undercool',
    'PDI': 1.45, 'npzfile':lala.npz}
    """
    with open(filename, 'r') as document:
        all_files_info_list = []
        pattern_chains_info = re.compile("c(\d+)n(\d+)")
        lines = document.readlines()
        lines = [line.strip() for line in lines]
        for line in lines:
            file_info_dict = {}
            chains_info_list = []
            print "the parse is true", line
        # if you need to parse through chain info
            if line.strip():  # non-empty line?
                # None means 'all whitespace', the default
                for match in pattern_chains_info.finditer(line):
                    chains_info_list.append((int(match.group(1)),
                                             int(match.group(2))))
            # tmp = re.findall("PDI(\d+).(\d+)", "lllPDI13.4")
            # tmp1 = re.findall(r"[-+]?\d*\.\d+|\d+", tmp[0])
            tmp = re.findall("PDI(\d+.\d+)extra", line)
            PDI = float(tmp[0])
            # PDI = 1.1
            npzfile = line
            extrastring = re.findall("extra(\S+).npz", line)
            file_info_dict['chains_info'] = chains_info_list
            file_info_dict['extrastring'] = extrastring[0]
            file_info_dict['PDI'] = PDI
            file_info_dict['npzfile'] = npzfile
            all_files_info_list.append(file_info_dict)
    print termcolor.colored(all_files_info_list, 'green')
    return all_files_info_list


def simple_parse_file_info(filename):
    """
    returns all_files_info_list with only basic information about files
    """
    with open(filename, 'r') as document:
        all_files_info_list = []
        lines = document.readlines()
        lines = [line.strip() for line in lines]
        for line in lines:
            file_info_dict = {}
            npzfile = line
            file_info_dict['npzfile'] = npzfile
            all_files_info_list.append(file_info_dict)
    print termcolor.colored(all_files_info_list, 'green')
    return all_files_info_list


def plot_points(args):
    """
    plots stuff with points
    #algorithm for induction time.
    filenames will represent the mol weight, and will be the X axis
    the flag = False will represent that the crystallization has started
    the tcryst = will be the time of crystallization, it will be appended
    to tcryst_array list
    """


    def plot_settings(parameters):
        plt.ylim(parameters['ymin'], parameters['ymax'])
        if parameters['xlog']:
            plt.xscale('log')
            parameters['xlabel'] = ''.join(('log', 'parameters[\'xlabel\']'))
            print "xlogscale"
        else:
            print "regular xscale"
        if parameters['ylog']:
            plt.yscale('log')
            parameters['ylabel'] = ''.join(('log', 'parameters[\'ylabel\']'))
            print "ylogscale"
        else:
            print "regular yscale"
        plt.ylabel('${0[ylabel]}$'.format(parameters))
        plt.xlabel('${0[xlabel]}$'.format(parameters))
        # plt.title('{0[title]}'.format(parameters))

        # if parameters['plot_style'] is 'lines':
        # points = ["-", "--", "-.", ":"]
        # elif parameters['plot_style'] is 'points':
        # points = [u'D', u'o', u'^', u'>', u's', u'8',
        #           u'<', u'>',  u'*',
        #           u'H', u'h', u'p', u'v', u'D',
        #           u'd', "-", "--", "-.", ":"]
        # # # else:
        points = [u'D-', u'o-', u'^-', u'>-', u's-', u'8-',
                  u'<-', u'>-',  u'*-',
                  u'H', u'h', u'p', u'v', u'D',
                  u'd', "-", "--", "-.", ":"]
        colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']
        lines = ["-", "--", "-.", ":"]
        return points, lines, colors

    # read file.param for the plot parameters and put them to parameters
    parameters = read_plot_parameters(args)

    # read file.in, maybe parse the chain info, reutrn the chain_info_parse
    chain_info_parse = parameters['parse_chain_info']
    if chain_info_parse is True:
        all_files_info_list = parse_file_info(args.input)
    elif chain_info_parse is False:
        all_files_info_list = simple_parse_file_info(args.input)
    else:
        raise ValueError('specify whether you need to parse chain info?')

    # initialize figure, do the plot parameters, and specify lines, ..
    plt.figure(figsize=(8, 7))
    points, lines, colors = plot_settings(parameters)
    linecycler = itertools.cycle(lines)
    pointscycler = itertools.cycle(points)
    colorcycler = itertools.cycle(colors)
    chi0_array, kappa_array, n_avrami_array, mw_array = [], [], [], []

    # loop through files listed and do the plotting
    for i, file_i in enumerate(all_files_info_list):
        # print "doing", file_i
        filename_i = all_files_info_list[i]['npzfile']

        # Extra_string = int(re.search('ICC(\s)(\d+)n(\d+)', j).group(1))
        npz = np.load(filename_i)
        t = npz['arr_0']
        p = npz['arr_1']
        if parameters['scale_dump'] is True:
            # then scale the result data
            # so the final data represents to the 10^6 lj untis values
            t *= parameters['timestep']*parameters['dumpskip']/1e6
        else:
            print "Using raw data without time scaling"
        p = p - p[0]
        # print p
        print t

        clr, llr, plr = next(colorcycler), next(linecycler), next(pointscycler)

        if parameters['parse_chain_info'] is True:
            fileNinfo_i = all_files_info_list[i]['chains_info']
            extrainfo_i = all_files_info_list[i]['extrastring']
            # pdiinfo_i = all_files_info_list[i]['PDI']
            result_name = ''
            for Mn, Nch in fileNinfo_i:
                result_name += 'M_n = '
                result_name += str(Mn)
                result_name += 'N_{chains} = '
                result_name += str(Nch)
                result_name += ';'
            # label_i = '${result_name}\ PDI={PDI};{extra}$'.format(
                # result_name=result_name, PDI=pdiinfo_i, extra=extrainfo_i)
            label_i = '${result_name}{extra}$'.format(
                result_name=result_name, extra=extrainfo_i)
            # label_i = '${result_name}$'.format(
            #     result_name=result_name)
        else:
            result_name = all_files_info_list[i]['npzfile']
            result_name = get_path_names.get_filename(result_name)
            label_i = '${result_name}$'.format(
                result_name=result_name)

        if parameters['fit'] is True:
            #  fitting
            #  deltaf/f <= sqrt(deltachi/chi^2 + chi0^2deltakappa/kappa
            #  chi0^2deltan/n)
            #  chi0 < = 1, thererfore i need just module of the whole thing
            #  this is the estimate of the partial derivatives df/dx
            #
            print "fitting data"
            if parameters['parse_chain_info'] is False:
                raise ValueError('specify that\
                 you need to parse chain info, to get M_w information')
            opt_p, pcov_p = scipy.optimize.curve_fit(func, t, p)
            # calculate error of each parameter using variance == diagonal
            # elements of covariance matrix
            chi0_err, kappa_err, n_avrami_err = np.sqrt(np.diag(pcov_p))
            tot_fit_err = np.linalg.norm(np.diag(pcov_p)/opt_p)
            print "total error %f" % tot_fit_err
            print "errors", np.sqrt(np.diag(pcov_p))
            chi0, kappa, n_avrami = opt_p
            chi0_array.append(chi0)
            kappa_array.append(kappa)
            n_avrami_array.append(n_avrami)
            mw, _ = fileNinfo_i[0]  # append M_w info to the list
            mw_array.append(mw)
            print "saturation", chi0
            print "kappa ", kappa

        # if parameters['normalize'] is True:
        #     print "normalizing data"
        #     # if parameters['fit'] is True:
        #     #     raise ValueError('fit and normalize are non-compatible')
        #     p /= p.max()
        #     if chain_info_parse is True:
        #         plt.plot(t, p, plr,
        #                  color=clr,
        #                  linewidth=1.8,
        #                  alpha=0.9,
        #                  label=label_i)
        #     else:
        #         plt.plot(t, p, llr,
        #                  color=clr,
        #                  linewidth=2.3,
        #                  alpha=0.95,
        #                  label=label_i)
        if parameters['fit'] is True:
            if parameters['normalize'] is True:
                raise ValueError('fit and normalize are non-compatible')
            plt.plot(t, func(t, chi0, kappa, n_avrami),
                     '-',
                     color=clr,
                     linewidth=4.2)
            plt.fill_between(t,
                             func(t, chi0, kappa, n_avrami)*(1.-tot_fit_err),
                             func(t, chi0, kappa, n_avrami)*(1.+tot_fit_err),
                             alpha=0.08, color=clr)
            if chain_info_parse is True:
                plt.plot(t, p, plr,
                         color=clr,
                         linewidth=1.8,
                         alpha=0.7,
                         label=label_i)
            else:
                plt.plot(t, p, llr,
                         color=clr,
                         linewidth=2.0,
                         alpha=0.95,
                         label=label_i)
        else:
            if parameters['normalize'] is True:
                print "normalizing data"
                p /= p.max()
            if chain_info_parse is True:
                plt.plot(t, p, plr,
                         color=clr,
                         linewidth=1.8,
                         alpha=0.6,
                         label=label_i)
            else:
                plt.plot(t, p, plr,
                         color=clr,
                         linewidth=2.8,
                         alpha=0.9,
                         label=label_i)
                # plt.plot(t, p, llr,
                #          color=clr,
                #          linewidth=2.0,
                #          alpha=0.95,
                #          label=label_i)
    if parameters['legend'] is True:
        # plt.legend(loc='best', prop={'size': 12})
        plt.legend(loc='best')
    print parameters

    if parameters['fit'] is True:
        chi0_array, kappa_array, n_avrami_array, mw_array =\
                                         np.asarray(chi0_array),\
                                         np.asarray(kappa_array),\
                                         np.asarray(n_avrami_array),\
                                         np.asarray(mw_array)
        print "chi0", chi0_array
        print "kappa", kappa_array
        print "n_avrami", n_avrami_array
        print "mw", mw_array
        np.savez('fit_data_chi0', mw_array, chi0_array)
        np.savez('fit_data_kappa', mw_array, kappa_array)
        np.savez('fit_data_n', mw_array, n_avrami_array)
        # np.savetxt('myfile.txt', np.c_[x,y,z])
        np.savetxt('parameters_{0[name]}.txt'.format(parameters),
                   np.c_[mw_array, chi0_array, kappa_array, n_avrami_array],
                   delimiter='\t', newline='\n', fmt='%.4f')
    plt.savefig('{0[plotname]}{0[name]}.pdf'.format(parameters))

    return None


def main():
    """main func"""
    args = read_parameters.read_from_file()
    # read_files('file.txt')
    # read_files(args)
    plot_points(args)
    return None

if __name__ == '__main__':
    main()

