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
import itertools

fig_width_pt = 400 # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width, fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': True,
          'figure.figsize': fig_size,
          "figure.autolayout": True}
plt.rcParams.update(params)

def fitting_function(params, t):
    """
    function: 1 - exp(-kappa*x**n),
    params: kappa, n (n_avrami)
    """
    # kappa, n = params
    kappa, t_0 = params
    # t_0 = t[np.where(x>0.01)[0][0]]
    # return (1. - np.exp(-kappa*(t-t_0)**3.))
    return (1. - np.exp(-kappa*(t-t_0)**3.))
    # return (1. - np.exp(-kappa*(t)**n))


def get_x0_extrapolation(data_t, data_x, extrapolation_name):
    """
    gets data of chi(t)
    this data will be extrapolated and chi(1/t) plotted
    intersection of which with 1/t = 0 will give the chi0
    """
    # print "data for exptrapolation", data_t, data_x
    inv_t = 1./data_t[(data_x>0.9*data_x.max())]
    f = np.polyfit(inv_t, data_x[(data_x>0.9*data_x.max())], 1)
    fit_fn = np.poly1d(f)
    # print "fit_fn = ", fit_fn
    fig_extrap = plt.figure()
    ax_extrap = fig_extrap.add_subplot(111)
    ax_extrap.set_xlabel('$1/(t-t_0)$')
    ax_extrap.set_ylabel('$\chi(1/(t-t_0))$')
    # t = np.linspace(0., inv_t.max(),0.1)
    t = np.linspace(0., inv_t.max() , 50)
    # ax_extrap.set_xlim([0, 15])
    # ax_extrap.plot(1./data_t[data_t>0.], data_x[data_t>0.], 'go',lw=2.)
    # # ax_extrap.plot(t, fit_fn(t),'r-')
    # ax_extrap.plot(t, fit_fn(t),'r-',lw=2.)
    # # ax_extrap.savefig('interpolation_x0.pdf')
    # fig_extrap.savefig(extrapolation_name + '.pdf')
    return fit_fn(0.)


# points = [u'D', u'o', u'^', u'>', u's', u'8',
#           u'<', u'>',  u'*',
#           u'H', u'h', u'p', u'v', u'D',
#           u'd', "-", "--", "-.", ":"]
points = ['D', 'o', '^', '>', 's', '8',
          '<', '>',  '*',
          'H', 'h', 'p', 'v', 'D',
          u'd', "-", "--", "-.", ":"]
colors = ['k', 'y', 'm', 'c', 'b', 'g', 'r', '#aaaaaa']
lines = ["-", "--", "-.", ":"]
linecycler = itertools.cycle(lines)
pointscycler = itertools.cycle(points)
colorcycler = itertools.cycle(colors)

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
chi0_array, chis_array, kappa_array, mw_array, n_avrami_array = [], [], [], [], []
chis_error_array = []

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
# loop through files listed and do the plotting
for i, file_i in enumerate(all_files_info_list):
    clr, llr, plr = next(colorcycler), next(linecycler), next(pointscycler)

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
    # t0 = t[x > 0.01][0]  # where the crystallization starts
    # xs = x.max()
    xs = get_x0_extrapolation(t, x, 'exptrapolation_xs')
    x0 = x.max()
    chi0_array.append(x0)
    chis_array.append(xs)
    chis_error_array.append(abs(xs-x0))
    mw_array.append(mw)
    x /= xs
    # t -= t0
    # print t
    # print x

    # # x /= x0
    ax2.set_xlabel(r'$\frac{t-t_0}{\tau}$')
    # ax2.set_ylabel('$\chi/\chi_s$')
    ax2.set_ylabel('$\chi(t)/\chi_s$')
    # ax2.set_xlabel('$t,\ 10^6\ lj$')
    # t0 = t[np.where(x>0.01)[0][0]]
    ax2.set_ylim([0, 1.1])
    ax2.set_xlim([0, 0.05])
    # ax2.set_xlim([0, 0.3])

    # print "t0 = ", t0
            # (tfitted, xfitted), fit_params, fit_params_err, fit_chi = fit.fit(fitting_function, t, x, default_pars = [20., 3.])
    (tfitted, xfitted), fit_params, fit_params_err, fit_chi = fit.fit(fitting_function, t[(x>0.01)&(x<0.6)], x[(x>0.01)&(x<0.6)], default_pars = [10.,0.1])

    kappa_array.append(fit_params[0])
    # n_avrami_array.append(fit_params[1])

    dummy_t = np.linspace(0.,0.6,100)
    # ax2.set_xscale('log')
    # ax2.set_xlim((0.01,0.3 ))
    if len(t) > 130:
        mrkevery = 2
    else:
        mrkevery = 1
    # ax2.plot(t, x, marker=plr,
    #     markersize=5,
    #      color=clr,
    #      alpha=0.95,
    #      label='n='+str(mw),
    #      markevery=mrkevery,
    #      zorder = 1,
    #      fillstyle=u'none')
    # # ,markevery=slice(0,10,2)
    # ax2.plot(dummy_t, fitting_function(fit_params, dummy_t ), '-',
    #      color=clr,
    #      linewidth=2.5,
    #      alpha=0.95,
    #      zorder = 10)
    # ax2.plot(t[x>0.01], x[x>0.01], 'o')
    # ax2.plot(dummy_t, fitting_function(fit_params, dummy_t ), ,lw=2.)
    ax2.plot((t[x>0.01]-fit_params[1])*fit_params[0]**(-1./3.), x[x>0.01],lw=2.,label='n='+str(mw))
ax2.legend(loc='best')
# fig2.savefig(filename_i+'fitting.pdf')
fig2.savefig('mastercurve_cool2slow.pdf')
#     chi0_array.append(chi0)
#     mw_array.append(mw)

chi0_array, mw_array =\
                                     np.asarray(chi0_array),\
                                     np.asarray(mw_array)
chis_array = np.asarray(chis_array)
chis_error_array = np.asarray(chis_error_array)
ind = np.argsort(mw_array)
print ind
print "mw", mw_array
# print "chi0", chi0_array
print "chis", chis_array
print "+-", chis_error_array
print "kappa", kappa_array
# print "n_avrami", n_avrami_array
# print "n_avrami", n_avrami_array
# print "mw", mw_array
# np.savez('x0'+extrainfo_i, mw_array[ind], chi0_array[ind])
np.savez('xs'+extrainfo_i, mw_array[ind], chis_array[ind],
 chis_error_array[ind])

#  to do :

#  output of chi(1/t) nice plot so people see how this shit is done
#  how do we determine max crystallinity value

# output of the values \kappa(Mw) -> tau(Mw)

#  master curve chi/chi_max (t/tau) which will comprise of every curve

# nice proper value plots chi(t) with fitting function provided

# output how was fitting done, so images of the fitted values so people see how bad or good was the fitting



