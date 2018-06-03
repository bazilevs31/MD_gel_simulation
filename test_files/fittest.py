#!/usr/bin/env python

import numpy as np
import itertools
import re
import math
# from lmfit import minimize
# from lmfit import minimize, Parameters, Parameter, report_errors,report_fit,conf_interval, printfuncs
import lmfit
# import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import termcolor
import read_parameters
import PlotFromFile


def residual(params, x, data=None):
    """
    model:
    alpha = \chi/chi0
    y = ln[-ln(1-alpha)]
    x = ln(t)
    y = n_avrami*x + lnkappa
    """
    lnkappa = params['lnkappa'].value             # parameter
    n_avrami = params['n_avrami'].value             # parameter
    # b = params['b'].value             # parameter
    model = lnkappa + x*n_avrami    # 2D function
    if data is None:
        return model
    return model - data


def fit_all(t, p):
    params = lmfit.Parameters()
    params.add('lnkappa', value=2.)
    params.add('n_avrami', value=1.)
    out = lmfit.minimize(residual, params, args=(t, p, ))  # lmfit minimizer
    lmfit.report_fit(params)

    # return out.params['chi0'].value
    return out


def renormalize_arrays(t, p):
    """returns x,y from t,p"""
    p = p - p[0]
    p0 = p.max()
    p = p/p0
    # print t, p
    #changed here
    # fitted_indices = np.where(np.logical_and(p >= 0.52, p <= 0.99))
    fitted_indices = np.where(np.logical_and(p >= 0.2, p <= 0.7))
    print fitted_indices
    p = p[fitted_indices]
    t = t[fitted_indices]
    # new coordinates
    x, y = np.log(t), np.log(-np.log(1-p))
    return x, y


def plot_settings(parameters, ax_plot):
    ax_plot.set_ylim(parameters['ymin'], parameters['ymax'])
    if parameters['xlog']:
        ax_plot.set_xscale('log')
        parameters['xlabel'] = ''.join(('log', 'parameters[\'xlabel\']'))
        print "xlogscale"
    else:
        print "regular xscale"
    if parameters['ylog']:
        ax_plot.set_yscale('log')
        parameters['ylabel'] = ''.join(('log', 'parameters[\'ylabel\']'))
        print "ylogscale"
    else:
        print "regular yscale"
    ax_plot.set_ylabel('${0[ylabel]}$'.format(parameters))
    ax_plot.set_xlabel('${0[xlabel]}$'.format(parameters))
    ax_plot.set_title('{0[title]}'.format(parameters))

    # if parameters['plot_style'] is 'lines':
    # points = ["-", "--", "-.", ":"]
    # elif parameters['plot_style'] is 'points':
    points = [u'D', u'o', u'^', u'>', u's', u'8',
              u'<', u'>',  u'*',
              u'H', u'h', u'p', u'v', u'D',
              u'd', "-", "--", "-.", ":"]
    # # else:
    # points = [u'D-', u'o-', u'^-', u'>-', u's-', u'8-',
    #           u'<', u'>',  u'*',
    #           u'H', u'h', u'p', u'v', u'D',
    #           u'd', "-", "--", "-.", ":"]
    colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']
    lines = ["-", "--", "-.", ":"]
    return points, lines, colors, ax_plot


def get_errors(fitted_parameters):
    """
    get errors by parsing a string of the output
    """
    printing_info = str(fitted_parameters.params['lnkappa'])
    tmp = re.findall("\s(\d+.\d+),", printing_info)
    kappa_error = float(tmp[0])
    printing_info = str(fitted_parameters.params['lnkappa'])
    tmp = re.findall("\s(\d+.\d+),", printing_info)
    n_avrami_error = float(tmp[0])
    return np.exp(kappa_error), n_avrami_error


def main():
    """
    main function
    Avrami Fitting:\\
    $\chi = \chi_0 (1 - e^{-\kappa t^n})$, \\
    $\alpha = \chi / \chi_0$,\\
    New Coordinates \\
    $x = ln(t)$, where $t: \alpha > 0.5$.\\
    $y = ln (-ln(1-\alpha))$\\
    So:\\
    $ y = ln(\kappa) + n x$
    """

    args = read_parameters.read_from_file()
    # read_files('file.txt')
    # read_files(args)
    # read file.param for the plot parameters and put them to parameters
    parameters = PlotFromFile.read_plot_parameters(args)

    # read file.in, maybe parse the chain info, reutrn the chain_info_parse
    chain_info_parse = parameters['parse_chain_info']
    if chain_info_parse is True:
        all_files_info_list = PlotFromFile.parse_file_info(args.input)
    elif chain_info_parse is False:
        all_files_info_list = PlotFromFile.simple_parse_file_info(args.input)
    else:
        raise ValueError('specify whether you need to parse chain info?')

    # initialize figure, do the plot parameters, and specify lines, ..
    chi0_array, kappa_array, n_avrami_array, mw_array = [], [], [], []
    kappa_error_array, n_avrami_error_array = [], []

    fig_plot = plt.figure(figsize=(8, 7))
    ax_plot = fig_plot.add_subplot(111)
    points, lines, colors, ax_plot = plot_settings(parameters, ax_plot)
    linecycler = itertools.cycle(lines)
    pointscycler = itertools.cycle(points)
    colorcycler = itertools.cycle(colors)

    # loop through files listed and do the plotting
    for i, file_i in enumerate(all_files_info_list):
        clr, llr, plr = next(colorcycler), next(linecycler), next(pointscycler)
        # print "doing", file_i
        filename_i = all_files_info_list[i]['npzfile']
        fileNinfo_i = all_files_info_list[i]['chains_info']
        extrainfo_i = all_files_info_list[i]['extrastring']
        result_name = ''
        result_name = ''
        for Mn, Nch in fileNinfo_i:
            result_name += 'M_n = '
            result_name += str(Mn)
            result_name += 'N_{chains} = '
            result_name += str(Nch)
            result_name += ';'
        # label_i = '${result_name}\ PDI={PDI};{extra}$'.format(
            # result_name=result_name, PDI=pdiinfo_i, extra=extrainfo_i)
        label_i = '${result_name}{extra}$'.format(result_name=result_name,
                                                  extra=extrainfo_i)
        mw, _ = fileNinfo_i[0]
        # Extra_string = int(re.search('ICC(\s)(\d+)n(\d+)', j).group(1))
        npz = np.load(filename_i)
        # if chain_info_parse is True:
        # else:
        #     ax_plot.plot(t, p, llr,
        #              color=clr,
        #              linewidth=2.3,
        #              alpha=0.95,
        #              label=label_i)
            # ax_plot.plot(t, p, 'go')
        # ax_plot.savefig('./figures/'+str(i)+'.pdf')
        # t = npz['arr_0']
        # p = npz['arr_1']
        # npz = np.load('testquench.npz')
        # npz = np.load('test.npz')
        t, p = npz['arr_0'], npz['arr_1']
        x, y = renormalize_arrays(t, p)
        chi0 = p.max()
        out = fit_all(x, y)
        kappa_error, n_avrami_error = get_errors(out)
        kappa_error_array.append(kappa_error)
        n_avrami_error_array.append(n_avrami_error)
        lnkappa, n_avrami = out.params['lnkappa'].value,\
                            out.params['n_avrami'].value
        tot_fit_err = np.sqrt((kappa_error/np.exp(lnkappa))**2 + \
                      (n_avrami_error/n_avrami)**2)
        ax_plot.plot(x, y, plr,
                     color=clr,
                     linewidth=1.6,
                     alpha=0.8,
                     label=label_i)
        ax_plot.fill_between(x,
                 (lnkappa + n_avrami*x)*(1.-tot_fit_err),
                 (lnkappa + n_avrami*x)*(1.+tot_fit_err),
                 alpha=0.08, color=clr)
        ax_plot.plot(x, lnkappa + n_avrami*x,
                     color=clr,
                     linewidth=2.8,
                     alpha=1.)
        mw_array.append(mw)
        kappa_array.append(np.exp(lnkappa))
        n_avrami_array.append(n_avrami)
        chi0_array.append(chi0)

    ax_plot.legend(loc='best')
    fig_plot.savefig('./figures/{0[plotname]}{0[name]}.pdf'.format(parameters))


    mw_array, chi0_array = zip(*sorted(zip(mw_array, chi0_array)))
    mw_array, kappa_array, kappa_error_array = zip(*sorted(zip(mw_array, kappa_array, kappa_error_array)))
    mw_array, n_avrami_array, n_avrami_error_array= zip(*sorted(zip(mw_array, n_avrami_array, n_avrami_error_array)))
    chi0_array, kappa_array, n_avrami_array, mw_array =\
                                     np.asarray(chi0_array),\
                                     np.asarray(kappa_array),\
                                     np.asarray(n_avrami_array),\
                                     np.asarray(mw_array)
    print "chi0", chi0_array
    print "kappa", kappa_array
    print "n_avrami", n_avrami_array
    print "mw", mw_array
    print "n_avrami_error_array", n_avrami_error_array

    # np.savez('./figures/fit{extra}chi0'.format(extra=extrainfo_i), mw_array, chi0_array)
    # np.savez('./figures/fit{extra}kappa'.format(extra=extrainfo_i), mw_array, kappa_array,kappa_error_array)
    # np.savez('./figures/fit{extra}n_avrami'.format(extra=extrainfo_i), mw_array, n_avrami_array, kappa_error_array)

    np.savez('./figures/fit{0[name]}chi0'.format(parameters), mw_array, chi0_array)
    np.savez('./figures/fit{0[name]}kappa'.format(parameters), mw_array, kappa_array,kappa_error_array)
    np.savez('./figures/fit{0[name]}n_avrami'.format(parameters), mw_array, n_avrami_array, n_avrami_error_array)

if __name__ == '__main__':
    main()
