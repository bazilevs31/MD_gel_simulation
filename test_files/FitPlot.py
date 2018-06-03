#!/usr/bin/env python

import numpy as np
# from lmfit import minimize
from lmfit import minimize, Parameters, Parameter, report_errors,report_fit,conf_interval, printfuncs
import matplotlib.pyplot as plt
import os
import re
# import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import termcolor
import read_parameters
# import scipy.optimize
import ast


def residual(params, x, data=None):
    chi0 = params['chi0'].value             # parameter
    kappa = params['kappa'].value             # parameter
    n_avrami = params['n_avrami'].value             # parameter
    # b = params['b'].value             # parameter
    model = chi0 * (1. - np.exp(-kappa * x**n_avrami))     # 2D function
    if data is None:
        return model
    return model - data


def fit_all(t, p, i):
    params = Parameters()
    params.add('chi0', value=0.5)
    params.add('kappa', value=2.)
    params.add('n_avrami', value=1.)
    out = minimize(residual, params, args=(t, p, ))  # lmfit minimizer
    report_fit(params)
    plt.plot(t, p, 'go')
    plt.plot(t, residual(out.params, t), 'r-')
    plt.savefig('./figures/'+str(i)+'.pdf')
    return out.params['chi0'].value

args = read_parameters.read_from_file()
# read_files('file.txt')
# read_files(args)
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
chi0_array, kappa_array, n_avrami_array, mw_array = [], [], [], []

# loop through files listed and do the plotting
for i, file_i in enumerate(all_files_info_list):
    # print "doing", file_i
    filename_i = all_files_info_list[i]['npzfile']
    fileNinfo_i = all_files_info_list[i]['chains_info']
    extrainfo_i = all_files_info_list[i]['extrastring']
    result_name = ''
    # print filename_i
    mw, _ = fileNinfo_i[0]

    # Extra_string = int(re.search('ICC(\s)(\d+)n(\d+)', j).group(1))
    npz = np.load(filename_i)
    t = npz['arr_0']
    p = npz['arr_1']
    p = p - p[0]
    chi0 = fit_all(t, p, i)
    chi0_array.append(chi0)
    mw_array.append(mw)

mw_array, chi0_array = zip(*sorted(zip(mw_array,chi0_array)))
chi0_array, mw_array =\
                                     np.asarray(chi0_array),\
                                     np.asarray(mw_array)
print "chi0", chi0_array
print "kappa", kappa_array
print "n_avrami", n_avrami_array
print "mw", mw_array
np.savez(extrainfo_i, mw_array, chi0_array)


