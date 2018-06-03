#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import plot_read_file
# import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import termcolor
import read_parameters
# import scipy.optimize
import ast
import fit


def get_x0_extrapolation(data_t, data_x, extrapolation_name):
    """
    gets data of chi(t)
    this data will be extrapolated and chi(1/t) plotted
    intersection of which with 1/t = 0 will give the chi0
    """
    print "data for exptrapolation", data_t, data_x
    inv_t = 1./data_t[(data_x>0.9*data_x.max())]
    f = np.polyfit(inv_t, data_x[(data_x>0.9*data_x.max())], 1)
    fit_fn = np.poly1d(f)
    print "fit_fn = ", fit_fn
    fig_extrap = plt.figure()
    ax_extrap = fig_extrap.add_subplot(111)
    # t = np.linspace(0., inv_t.max(),0.1)
    t = np.linspace(0., inv_t.max() , 50)
    ax_extrap.set_xlim([0, 5])
    ax_extrap.plot(1./data_t[data_t>0.], data_x[data_t>0.], 'go')
    # ax_extrap.plot(t, fit_fn(t),'r-')
    ax_extrap.plot(t, fit_fn(t),'r-')
    # ax_extrap.savefig('interpolation_x0.pdf')
    fig_extrap.savefig(extrapolation_name + '.pdf')
    return fit_fn(0.)



args = read_parameters.read_from_file()
parameters = plot_read_file.read_plot_parameters(args)
chain_info_parse = parameters['parse_chain_info']

if chain_info_parse is True:
    all_files_info_list = plot_read_file.parse_file_info(args.input)
elif chain_info_parse is False:
    all_files_info_list = plot_read_file.simple_parse_file_info(args.input)
else:
    raise ValueError('specify whether you need to parse chain info?')

# initialize figure, do the plot parameters, and specify lines, ..
# chi0_array, kappa_array, mw_array = [], [], [], []


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
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
    x = npz['arr_1']
    x = x - x[0]
    t0 = t[x > 0.01][0]  # where the crystallization starts
    # xs = x.max()
    xs = get_x0_extrapolation(t, x, 'exptrapolation_xs')
    x /= xs
    # t -= t0
    print t
    print x
    new_x = np.log(-np.log(1.-x[(x>0.01)&(x<0.3)]))
    new_t = np.log(t[(x>0.01)&(x<0.3)])
    print new_t
    print new_x
    f = np.polyfit(new_t, new_x, 1)
    fit_fn = np.poly1d(f)
    print fit_fn
    ax2.plot(new_t, new_x, 'go', new_t, fit_fn(new_t), 'k-')
fig2.savefig('fitting.pdf')
#     chi0_array.append(chi0)
#     mw_array.append(mw)

# mw_array, chi0_array = zip(*sorted(zip(mw_array,chi0_array)))
# chi0_array, mw_array =\
#                                      np.asarray(chi0_array),\
#                                      np.asarray(mw_array)
# print "chi0", chi0_array
# print "kappa", kappa_array
# print "n_avrami", n_avrami_array
# print "mw", mw_array
# np.savez(extrainfo_i, mw_array, chi0_array)

#  to do :

#  output of chi(1/t) nice plot so people see how this shit is done
#  how do we determine max crystallinity value

# output of the values \kappa(Mw) -> tau(Mw)

#  master curve chi/chi_max (t/tau) which will comprise of every curve

# nice proper value plots chi(t) with fitting function provided

# output how was fitting done, so images of the fitted values so people see how bad or good was the fitting



