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
import itertools
import numpy as np

import termcolor
import read_parameters
import get_path_names
import scipy.optimize
import plot_read_file
import matplotlib.pyplot as plt
import fit


def get_x0_extrapolation(data_t, data_x):
    """
    gets data of chi(t)
    this data will be extrapolated and chi(1/t) plotted
    intersection of which with 1/t = 0 will give the chi0
    """
    inv_t = np.power(data_t, -1.)
    fit = np.polyfit(inv_t, data_x, 1)
    fit_fn = np.poly1d(fit)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(1./inv_t, data_x, 'go',1./inv_t, fit_fn(1./inv_t),'k-')
    # ax1.savefig('interpolation_x0.pdf')
    plt.savefig('interpolation_x0.pdf')
    return fit_fn(0.)

def fitting_function(params, t):
    """
    function: 1 - exp(-kappa*x**n),
    params: kappa, n (n_avrami)
    """
    # kappa, n = params
    kappa = params
    return (1. - np.exp(-kappa*(t)**3.))
    # return (1. - np.exp(-kappa*(t)**n))



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
    parameters = plot_read_file.read_plot_parameters(args)

    # read file.in, maybe parse the chain info, reutrn the chain_info_parse
    chain_info_parse = parameters['parse_chain_info']
    if chain_info_parse is True:
        all_files_info_list = plot_read_file.parse_file_info(args.input)
    elif chain_info_parse is False:
        all_files_info_list = plot_read_file.simple_parse_file_info(args.input)
    else:
        raise ValueError('specify whether you need to parse chain info?')

    # initialize figure, do the plot parameters, and specify lines, ..
    plt.figure(figsize=(8, 7))
    points, lines, colors = plot_settings(parameters)
    linecycler = itertools.cycle(lines)
    pointscycler = itertools.cycle(points)
    colorcycler = itertools.cycle(colors)
    chi0_array, kappa_array, n_avrami_array, mw_array = [], [], [], []
    chi0_array_error = []
    chis_array = []
    # loop through files listed and do the plotting
    for i, file_i in enumerate(all_files_info_list):
        # print "doing", file_i
        filename_i = all_files_info_list[i]['npzfile']

        npz = np.load(filename_i)
        t = npz['arr_0']
        x = npz['arr_1']
        x -= x[0]
        # t =
        # normalizing t to be from 0, 0.3
        # t /= t.max()
        # t *= 0.3
        print x
        if parameters['scale_dump'] is True:
            # then scale the result data
            # so the final data represents to the 10^6 lj untis values
            t *= parameters['timestep']*parameters['dumpskip']/1e6
        else:
            print "Using raw data without time scaling"
        #x = x - x[0]
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
        mw, _ = fileNinfo_i[0]  # append M_w info to the list
        mw_array.append(mw)
        x0 = x.max()
        xs = max(get_x0_extrapolation(t[-10:], x[-10:]),x.max())
        chi0_array.append(x0)
        chis_array.append(xs)
        chi0_array_error.append(abs(xs-x0))

        # if parameters['fit'] is True:
        #     #  fitting
        #     #  deltaf/f <= sqrt(deltachi/chi^2 + chi0^2deltakappa/kappa
        #     #  chi0^2deltan/n)
        #     #  chi0 < = 1, thererfore i need just module of the whole thing
        #     #  this is the estimate of the partial derivatives df/dx
        #     #
        #     print "fitting data"
        #     if parameters['parse_chain_info'] is False:
        #         raise ValueError('specify that\
        #          you need to parse chain info, to get M_w information')
        #
        #     # x0 = x.max()
        #     x0 = max(get_x0_extrapolation(t[-10:], x[-10:]),x.max())
        #     # plt.clf()
        #     x0_err = 0.
        #     x /= x0
        #     t0 = t[np.where(x>0.01)[0][0]]
        #     print "t0 = ", t0
        #     # (tfitted, xfitted), fit_params, fit_params_err, fit_chi = fit.fit(fitting_function, t, x, default_pars = [20., 3.])
        #     (tfitted, xfitted), fit_params, fit_params_err, fit_chi = fit.fit(fitting_function, t[x>0.01], x[x>0.01], default_pars = [10.])
        #
        #
        #     print "total error %f" % fit_chi
        #     # print "errors", np.sqrt(np.diag(pcov_p))
        #     print "x0 = ", x0, " +- ", x0_err
        #     print "kappa = ", fit_params[0], " +- ", fit_params_err[0]
        #     # print "n_avrami = ", fit_params[1], " +- ", fit_params_err[1]
        #     chi0_array.append(x0)
        #     kappa_array.append(fit_params[0])
        #     # n_avrami_array.append(fit_params[1])
        #     mw, _ = fileNinfo_i[0]  # append M_w info to the list
        #     mw_array.append(mw)

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


        if chain_info_parse is True:
            plt.plot(t, x, llr,
                     color=clr,
                     linewidth=2.4,markevery=100,
                     alpha=0.99,
                     label="n = " + str(mw))

        if parameters['fit'] is True:
            if parameters['normalize'] is True:
                raise ValueError('fit and normalize are non-compatible')
            plt.plot(t, fitting_function(fit_params, t),
                     '-',
                     color=clr,
                     linewidth=4.2)
            # plt.fill_between(t,
            #                  func(t, chi0, kappa, n_avrami)*(1.-tot_fit_err),
            #                  func(t, chi0, kappa, n_avrami)*(1.+tot_fit_err),
            #                  alpha=0.08, color=clr)
            if chain_info_parse is True:
                plt.plot(t, x, plr,
                         color=clr,
                         linewidth=1.8,
                         alpha=0.7,
                         label=label_i)
            else:
                plt.plot(t, x, llr,
                         color=clr,
                         linewidth=2.0,
                         alpha=0.95,
                         label=label_i)
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
        # np.savetxt('parameters_{0[name]}.txt'.format(parameters),
                   # np.c_[mw_array, chi0_array, kappa_array, n_avrami_array],
                   # delimiter='\t', newline='\n', fmt='%.4f')
    chi0_array, chis_array, chi0_array_error, mw_array =\
                                         np.asarray(chi0_array),\
                                         np.asarray(chis_array),\
                                         np.asarray(chi0_array_error),\
                                         np.asarray(mw_array)

    np.savetxt(extrainfo_i+'fit_data_Mn_chi0_delta.csv', np.c_[mw_array, chi0_array,chis_array, chi0_array_error],delimiter=',')

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

