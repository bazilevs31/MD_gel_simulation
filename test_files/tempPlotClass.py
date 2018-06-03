
import os
import itertools
import re
import numpy as np
# from lmfit import minimize, Parameters, Parameter, report_errors,report_fit,conf_interval, printfuncs
import lmfit
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import termcolor
import read_parameters
import get_path_names
import scipy.optimize
import ast

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
                           'parse_chain_info']:
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



class PlotFiles(object):
    """this class plots the .npz files in the .in files"""
    def __init__(self, x, y):
        # self.fig_plot = plt.figure()
        self.fig_norm = plt.figure()
        # self.fig_fit = plt.figure()
        # self.ax_plot = self.fig_plot.add_subplot(111)
        self.ax_norm = self.fig_norm.add_subplot(111)
        # self.ax_fit = self.fig_fit.add_subplot(111)
        self.x = x
        self.y = y

    def fit_all(self):
        """fit data in self variables"""
        def _residual(params, x, data=None):
            """
            model function for fitting
            """
            chi0 = params['chi0'].value             # parameter
            kappa = params['kappa'].value             # parameter
            n_avrami = params['n_avrami'].value             # parameter
            # b = params['b'].value             # parameter
            model = chi0 * (1. - np.exp(-kappa * x**n_avrami))     # 2D function
            if data is None:
                return model
            return model - data
        params = lmfit.Parameters()
        params.add('chi0', value=0.5)
        params.add('kappa', value=2.)
        params.add('n_avrami', value=1.)
        t = self.x
        p = self.y
        out = lmfit.minimize(_residual, params, args=(t, p, ))
        lmfit.report_fit(params)
        # plt.plot(t, p, 'go')
        self.ax_fit.plot(t, _residual(out.params, t), 'r-')
        # self.fig_fit.savefig('./figures/'+'fit.pdf')
        # return out.params['chi0'].value
        # yield eself.fig_fit
    def normalize_by_max(self):
        """takes an array of elements to plot
        gets maximum
        normalizes by it
        so the crystallinity is from 0 - 100 %
        """
        t = self.x
        # p = self.y/self.y.max()
        p = self.y
        print t,p
        self.ax_norm.plot(t, p, 'r--')
        # self.fig_fit.savefig('./figures/'+'fit.pdf')
        # then do the plotting
        # return None
        # /self.fig_norm
        yield self.fig_norm



    # def plot_func(self):

    #     self.ax1.plot(self.x, self.y)
    #     self.ax2.plot(self.x, self.y)
    #     self.fig1.savefig('./figures/test1.pdf')
    #     self.fig2.savefig('./figures/test2.pdf')

# p.plot_func()
x = np.arange(0.0, 2.0, 0.01)
for i in range(4):
    y = 0.6 * (1. - np.exp(-1.0 * x**i))
    p = PlotFiles(x, y)
    gen = p.normalize_by_max()
fig = p.fig_norm
fig.savefig('./figures/cool.pdf')
