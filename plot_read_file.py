#!/usr/bin/env python

# Program: plot_read_file.py
# Purpose: provides various functions for plotall type of programs
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python PlotFromFile.py for help,


# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


import os
import re
import termcolor
import ast


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

