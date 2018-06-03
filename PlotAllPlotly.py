#!/usr/bin/env python
import plotly.plotly as py
py.sign_in('vasiliy.triandafilidi', 'h5uejpbuf4')
import plotly.tools as tls

# (*) Graph objects to piece together plots
from plotly.graph_objs import *

import numpy as np  # (*) numpy for math functions and arrays
import matplotlib.pyplot as plt

import read_parameters
import get_path_names
import os
import itertools
import re
import plot_read_file

def plot_mpl_fig():
    """
    plots stuff with points
    #algorithm for induction time.
    filenames will represent the mol weight, and will be the X axis
    the flag = False will represent that the crystallization has started
    the tcryst = will be the time of crystallization, it will be appended
    to tcryst_array list
    """
    args = read_parameters.read_from_file()
    def plot_settings(parameters):
        # read_files('file.txt')
        # read_files(args)
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

    # loop through files listed and do the plotting
    for i, file_i in enumerate(all_files_info_list):
        # print "doing", file_i
        filename_i = all_files_info_list[i]['npzfile']

        npz = np.load(filename_i)
        t = npz['arr_0']
        x = npz['arr_1']
        is_with_error_bars = True
        try:
            error_array = npz['arr_2']
        except Exception:
            print "doing Exception"
            is_with_error_bars = False
            pass
        # t =
        print x
        if parameters['scale_dump'] is True:
            # then scale the result data
            # so the final data represents to the 10^6 lj untis values
            t *= parameters['timestep']*parameters['dumpskip']/1e6
        else:
            print "Using raw data without time scaling"
        # print p
        print t

        clr, llr, plr = next(colorcycler), next(linecycler), next(pointscycler)

        result_name = all_files_info_list[i]['npzfile']
        result_name = get_path_names.get_filename(result_name)
        label_i = '${result_name}$'.format(
            result_name=result_name)

        if is_with_error_bars:
            # plt.errorbar(t, x,
            #      yerr=error_array,
            #      color=clr,
            #      linewidth=1.8,
            #      label=label_i)
            plt.plot(t, x,
                 color=clr,
                 linewidth=1.8,
                 label=label_i)
            plt.errorbar(t, error_array,
                 color=clr,
                 linewidth=1.8,
                 label='error'+label_i)
        else:
            x = x - x[0]
            plt.plot(t, x, plr,
                     color=clr,
                     linewidth=1.8,
                     alpha=0.7,
                     label=label_i)
    print parameters
    return None

# Package all mpl plotting commands inside one function
# def plot_mpl_fig():

#     args = read_parameters.read_plot()
#     print args
#     filenames = [get_path_names.get_filename(x) for x in args.collection]
#     parameters = {'title': args.titlename,
#                   'xlabel': args.xlabel,
#                   'ylabel': args.ylabel,
#                   'dumpskip': args.dumpskip}

#     if args.logplot:
#         plt.xscale('log')
#         parameters['xlabel'] = ''.join(('log', 'parameters[\'xlabel\']'))
#         print "xlogscale"
#     else:
#         print "regular xscale"

#     lines = [u'D-', u'o-', u'^-', u'>-', u's', u'8', u'<', u'>',  u'*',
#              u'H', u'h', u'p', u'v', u'D', u'd', "-", "--", "-.", ":"]
#     # widths = [4, 6, 8, 10, 12]
#     colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']

#     linecycler = itertools.cycle(lines)
#     colorcycler = itertools.cycle(colors)

#     plt.figure()
#     plt.ylabel('${0[ylabel]}$'.format(parameters))
#     plt.xlabel(r'${0[xlabel]}$'.format(parameters))
#     # plt.title('${0[title]}$'.format(parameters))
#     # plt.ylim(0,1)
#     # plt.xlim(10,210)
#     count = 0
#     for (i, j) in itertools.izip(args.collection, filenames):
#         # ChainLength = int(re.search('(\d+)n(\d+)', j).group(1))
#         # NumOfChains = int(re.search('(\d+)n(\d+)', j).group(2))
#         print j
#         count += 1
#         # Extra_string = (re.findall('fit(\S+)', j))
#         npz = np.load(i)
#         t = npz[args.name_x]
#         p = npz[args.name_y]
#         # p_err = npz[args.name_z]
#         p_err = t*0.
#         # p = p - p[0]
#         clr = next(colorcycler)
#         # namefile=get_path_names.get_filename(j)
#         # namefile = Extra_string[0]
#         namefile = j
#         plt.plot(t, p,
#                      color=clr,
#                      ls='o-',
#                      linewidth=1.4,
#                      label=str(j))
#         # plt.plot(t, p, 'go-',color=clr,linewidth=1.4,label=namefile)
#         # plt.plot(t, p, next(linecycler),
#         #          color=clr,
#         #          linewidth=1.2,
#         #          label=namefile)
#                  #,
#                  # label=r'$M_n = %d,\ N_{chains} = %d$'
#                  # % (ChainLength, NumOfChains))
#     # if args.legend:
#     # plt.legend(loc='best')
#     # plt.savefig('{0[title]}.pdf'.format(parameters))
#     return None
# Plot it!

#use the function described above to plot everything
plot_mpl_fig()
# get current figure, retrieve the information on the figure that has been just created
#alternative mpl_fig1 = plt.figure()
mpl_fig1 = plt.gcf()

#now we created a plotly object and we are printing the information about it
py_fig1 = tls.mpl_to_plotly(mpl_fig1, verbose=True)

# print information
# print(py_fig1.to_string())
# print py_fig1
for i in ['autosize', 'width', 'height']:
    print i, py_fig1['layout'][i]

py_fig1['layout'].pop('annotations', None)

# Add legend, place it at the top right corner of the plot
py_fig1['layout'].update(
    showlegend=True,
    legend=Legend(
        x=1.05,
        y=1
    )
)

plot_url = py.plot_mpl(mpl_fig1, strip_style=True,
             filename='s6_damped_oscillation-default-style')
