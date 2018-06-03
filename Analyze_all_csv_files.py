#!/usr/bin/env python

import numpy as np
import pandas as pd
import re, os, argparse, pickle, termcolor
import Gel_dump_analysis

import Analyze_file_info
import matplotlib.pyplot as plt
# MATPLOTLIB_STYLE_USE = 'seaborn-poster'
MATPLOTLIB_STYLE_USE = 'presentation-new'
# MATPLOTLIB_STYLE_USE = 'presentation'
# MATPLOTLIB_STYLE_USE = 'seaborn'
# MATPLOTLIB_STYLE_USE = 'paper-doublefig'
plt.style.use(MATPLOTLIB_STYLE_USE)

def plot_local(x, y, y_divide_array_1, y_divide_array_2, ax, plot_details, plot_pure, plot_line=True, **kwargs):
    """plots the x, y array dividing by arrays

    Args:
        x (array): x
        y (array): y
        y_divide_array_1 (array): to divide
        y_divide_array_2 (array): to divide
        ax (matplotlib ax): to plot on
        plot_details (int): to which degree to plot
        plot_pure (bool): should I plot the pure data?
        **kwargs: fillstyle, color, filename

    Returns:
        None: Description
    """
    fillstyle = kwargs['fillstyle']
    c = kwargs['color']
    myfile = kwargs['filename']
    if fillstyle is not 'none':
        mfc = c
    else:
        mfc = fillstyle

    fillstyle = 'none'
    mfc = 'none'
    # since the values are logarithmic we will just need to
    # add/subtract
    if plot_pure:
        ax.plot(x, y,'o',fillstyle='none', mec=c, markersize=10.,markeredgewidth=1.3, label='$f^{eff}$_'+myfile)
    if plot_line:
        ax.plot(np.linspace(0,1,10), np.linspace(0,1,10), 'k--')
    if plot_details == 2:
        ax.plot(x, y/(y_divide_array_1),'o',fillstyle=fillstyle, mec=c,mfc=mfc, markersize=10.,markeredgewidth=1.3, label='$f^{eff}$_'+myfile)
    if plot_details == 3:
        ax.plot(x, y/(y_divide_array_1 - y_divide_array_2),'o',fillstyle=fillstyle, mec=c,mfc=mfc, markersize=10.,markeredgewidth=1.6, label='$f^{eff}$_'+myfile)
    if plot_details == 4:
        ax.plot(x, y/(y_divide_array_1),'o',fillstyle=fillstyle, mec=c,mfc=mfc, markersize=10.,markeredgewidth=1.3, label='$f^{eff}$_'+myfile)
        ax.plot(x, y/(y_divide_array_1 - y_divide_array_2),'>',fillstyle=fillstyle, mec=c,mfc=mfc, markersize=10.,markeredgewidth=1.3, label='$f^{eff}$_'+myfile)
    return None

def plot_ax(x, y, ax, **kwargs):
    """plots x and y with params

    Args:
        x (TYPE): Description
        y (TYPE): Description
        ax (TYPE): Description
        **kwargs: Description

    Returns:
        TYPE: Description
    """

    ax.set(xlim=[kwargs['xmin'], kwargs['xmax']],ylim=[kwargs['ymin'],kwargs['ymax']], xlabel='{0[xlabel]}'.format(kwargs), ylabel='{0[ylabel]}'.format(kwargs))
    ax.plot(x, y, kwargs['style'])
    return None

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.

    Args:
        n (TYPE): Description
        name (str, optional): Description

    Returns:
        TYPE: Description
    '''
    return plt.cm.get_cmap(name, n)

def in1d(a, b, tol=1e-3):
    """looks for similarities in two arrays with a certain tolerance

    Args:
        a (array): Description
        b (array): Description
        tol (float, optional): Description

    Returns:
        TYPE: Description
    """
    a = np.unique(a)
    intervals = np.empty(2*a.size, float)
    intervals[::2] = a - tol
    intervals[1::2] = a + tol
    overlaps = intervals[:-1] >= intervals[1:]
    overlaps[1:] = overlaps[1:] | overlaps[:-1]
    keep = np.concatenate((~overlaps, [True]))
    intervals = intervals[keep]
    return np.searchsorted(intervals, b, side='right') & 1 == 1

def get_arrays_from_file(**kwargs):
    """analyzes the file and returns the arrays with the arrays
    Args:
        **kwargs: Description

        temp_name_array
        ax_dict
        fig_dict
        bound_prefix
        bound_suffix
        myfile
        cond_prefix
        cond_suffix

    Returns:
        TYPE: Description
    """
    myfile = '{0[myfile]}'.format(kwargs)
    f_divide_array = []
    conc_divide_array = []
    phi3_array_pure = []
    temp_array = []
    press_array = []
    work_array = []
    pzz_array = []
    phi_unit_convert = 567.

    p_scale_factor = 1.

    # phi_press_mod_array = []
    # temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']

    # for f in files:
    temp_name_array = kwargs['temp_name_array']
    ax_dict = kwargs['ax_dict']
    fig_dict = kwargs['fig_dict']

    for temp_name in temp_name_array:
        f = myfile+'_t'+temp_name
        if os.path.exists(f+'.csv') and ('fl01_fr02_t01' not in f):
            fl, fr, temp, gelname = Analyze_file_info.get_profile_params(f)
            temp_array.append(temp)
            f += '.csv'
            df = pd.read_csv(f)

            phi_max, phi_min = Analyze_file_info.get_steps_from_profile(df['coord_phi'],df['phiz'], ax_dict['phiz'], fig_dict['phiz'], variable='phiz', name=f)

            phi3_array_pure.append(np.abs(phi_max - phi_min)/(phi_unit_convert))

            if kwargs['press_normalize_ncount'] is False:
                press_max, press_min = Analyze_file_info.get_steps_from_profile(df['Coord1'],df['Pxxyy'], ax_dict['Plat'], fig_dict['Plat'], variable='Plat', name=f)
            elif kwargs['press_normalize_ncount'] is True:
                press_max, press_min = Analyze_file_info.get_steps_from_profile(df['Coord1'],df['Pxxyy']/df['Ncount'], ax_dict['Plat'], fig_dict['Plat'], variable='Plat', name=f)
            else:
                raise NameError('there was no press_normalize_ncount')

            press_array.append(p_scale_factor*np.abs(press_max - press_min))
            # press_array.append((press_r-press_l))

            concentration_max, concentration_min = Analyze_file_info.get_steps_from_profile(df['Coord1'],df['NumDensity_ci'], ax_dict['concentration'], fig_dict['concentration'], variable='concentration', name=f)
            conc_divide_array.append(np.log(concentration_max/concentration_min))


            pzz = np.mean(df['Pzz'][2:-2])
            pzz_array.append(pzz)

            f_divide_array.append(np.log(fr/fl))


    boundaryfile = '{0[bound_prefix]}_{0[myfile]}_{0[bound_suffix]}.csv'.format(kwargs)
    # boundaryfile = 'Boundary_' + myfile + '_interpolated' + '.csv'
    if os.path.exists(boundaryfile):
        print("Boundary file exists")
        print(boundaryfile)
        df_bound = pd.read_csv(boundaryfile)
        dl_array = df_bound['dl']
        V_array = 200*200*dl_array
        _t1_array = np.array(temp_array)
        _t2_array = np.linspace(0.1, 1.1, 11)
        # indices = np.where(np.in1d(_t2_array, _t1_array))[0]
        indices = np.nonzero(in1d(_t2_array, _t1_array))[0]

        print('indices')
        print(indices.shape)
        print(indices)
        print('pzz_array')
        print(np.array(pzz_array).shape)
        print('V_array')
        print(V_array.shape)
        print('temp_array')
        print(np.array(temp_array).shape)
        # print(temp_array)
        print(_t1_array)
        print(_t2_array)

        work_array =  np.array(pzz_array) * V_array[indices] # *

    # # if os.path.exists('COND_q_'+myfile+'_v2'+'.csv'):
    # csvfile = 'COND_q1_method_'+myfile+'_v2'+'.csv'
    # # csvfile = 'COND_v1_'+myfile+'_interpolated'+'.csv'

    csvfile = '{0[cond_prefix]}_{0[myfile]}_{0[cond_suffix]}.csv'.format(kwargs)
    if os.path.exists(csvfile):
        print("csvfile exists")
        print(csvfile)
        # df4 = pd.read_csv('COND_q_'+myfile+'.csv')
        df4 = pd.read_csv(csvfile)

        # print(temp_array)
        # divide_array = np.array(df4.loc[df4['Temp'].isin(iter(temp_array))])
        tt = np.array(temp_array)
        print('temp_array shape')
        print(tt.shape)
        tt2 = 10*tt-1
        tt3=tt2.astype(np.int32)
        divide_array = np.array(df4['param'])
        q1_q2_param_array = np.log(divide_array[tt3])

        fillstyle='full'
        # print(divide_array[tt3])
    else:
        # q1_q2_param_array = np.log(np.ones_like(np.array(temp_array)))
        q1_q2_param_array = np.log(np.ones_like(np.array(temp_array)))
        fillstyle='none'

    return np.array(temp_array), np.array(phi3_array_pure), np.array(conc_divide_array), np.array(f_divide_array), np.array(q1_q2_param_array), np.array(press_array), np.array(work_array),np.array(pzz_array), fillstyle

def main(args):
    """does the analysis

    Args:
        args (TYPE): Description
    """
    import matplotlib as mpl
    # mpl.rcParams.update(mpl.rcParamsDefault)
    list_of_input_axes = ['f','deltap','deltapT','c', 'work']
    list_of_intern_axes = ['phiz', 'Plat','concentration']
    list_of_axes = list_of_input_axes + list_of_intern_axes
    axes_index_list = range(1,len(list_of_axes)+1)

    index_dict = dict(zip(axes_index_list, list_of_axes))

    fig_dict = {}
    ax_dict = {}

    for fig_num, fig_str in index_dict.items():
        if fig_str in list_of_intern_axes:
            mpl.rcParams.update(mpl.rcParamsDefault)
            # plt.style.use('seaborn')
        elif fig_str in list_of_input_axes:
            plt.style.use(MATPLOTLIB_STYLE_USE)
        fig_dict[fig_str] = plt.figure(fig_num)
        ax_dict[fig_str]= fig_dict[fig_str].add_subplot(111)

    files_list = ['sc1844_nm100_l50_fl01_fr02','sc1844_nm100_l50_fl01_fr03','sc1844_nm100_l50_fl01_fr05','sc1844_nm100_l50_fl01_fr1','sc1844_nm100_l50_fl025_fr05']
    file_color_list = ['r','m','g','b','c']

    # files_list = ['sc1844_nm100_l50_fl01_fr02']
    # file_color_list = ['r','m','g','b','c']
    files_color_dict = dict(zip(files_list, file_color_list))

    temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
    temp_color_array = get_cmap(len(temp_name_array))

    params_dict = vars(args)

    params_dict['temp_name_array'] = temp_name_array
    params_dict['ax_dict'] = ax_dict
    params_dict['fig_dict'] = fig_dict
    # for myfile, c in files_color_dict.items():
    for myfile in params_dict['myfiles_list']:
        c = files_color_dict[myfile]
        params_dict['myfile'] = myfile
        temp_array, phi3_array_pure, conc_divide_array, f_divide_array, q1_q2_param_array, press_array,work_array, pzz_array, fillstyle = get_arrays_from_file(**params_dict)

        # Natoms = 101300
        # K_phi_P_array = (1/(2*Natoms*(conc_divide_array-q1_q2_param_array)))*((np.exp(conc_divide_array-q1_q2_param_array) - 1)/(np.exp(conc_divide_array-q1_q2_param_array )+1)

        for my_axes_ind, my_axes_name in index_dict.items():
            for my_param_name in params_dict.keys():
                if my_axes_name == my_param_name:
                    print('doing axis = {axis}'.format(axis=my_axes_name))
                    if params_dict[my_param_name]:
                        if my_param_name == 'f':
                            plot_local(temp_array, phi3_array_pure, f_divide_array, q1_q2_param_array, ax_dict[my_param_name], plot_details=params_dict['plot_details'],plot_pure=params_dict['plot_pure'],color=c,filename=myfile,fillstyle=fillstyle)
                        elif my_param_name == 'c':
                            plot_local(temp_array, phi3_array_pure, conc_divide_array, q1_q2_param_array, ax_dict[my_param_name], plot_details=params_dict['plot_details'],plot_pure=params_dict['plot_pure'],color=c,filename=myfile,fillstyle=fillstyle)
                        elif my_param_name == 'deltapT':
                            # plot_local(temp_array, press_array, np.ones_like(press_array), np.zeros_like(press_array), ax_dict[my_param_name], plot_details=1,plot_pure=True,color=c,filename=myfile,fillstyle=fillstyle)
                            plot_local(temp_array, press_array*200*200*1000, conc_divide_array, q1_q2_param_array, ax_dict[my_param_name], plot_details=params_dict['plot_details'],plot_pure=params_dict['plot_pure'],color=c,filename=myfile,fillstyle=fillstyle)
                        elif my_param_name == 'deltap':
                            p_scale_factor = 10**(int(params_dict['press_scale']))
                            # plot_local(press_array*200*200*1000/p_scale_factor, phi3_array_pure, conc_divide_array, np.zeros_like(press_array), ax_dict[my_param_name], plot_details=2,plot_pure=True, plot_line=False,color=c,filename=myfile,fillstyle=fillstyle)
                            plot_local(press_array*200*200*1000/p_scale_factor, phi3_array_pure, np.ones_like(phi3_array_pure), np.zeros_like(press_array), ax_dict[my_param_name], plot_details=1,plot_pure=True, plot_line=False,color=c,filename=myfile,fillstyle=fillstyle)


    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    # mathtext.fontset : cm
    # rc('text', usetex=True)
# MATPLOTLIB_STYLE_USE = 'presentation-new'
    plt.rcParams["font.family"] = "serif"
    # plt.rcParams["mathtext.fontset"] = "dejavuserif"
    plt.rcParams["mathtext.fontset"] = "stix"
    # plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams["font.serif"] = "STIXGeneral"
    for key, value in ax_dict.items():
        if key not in list_of_intern_axes:
        # fig_dict[key].show()
            # ax_dict[key].set(xlim=[0, 1.2], xlabel='$T \ [\\varepsilon]$', ylabel='$\\Delta \\varphi^{*} \ [\\varepsilon]$')
            ax_dict[key].set(xlabel=r'$T \ [\varepsilon]$', ylabel=r'$q\Delta \varphi^{*} \ [\varepsilon]$')

    # ax_dict['cond'].set(ylim=[0,3], ylabel='$ln \dfrac{(c_2 (1-q_2)} {c_1 (1-q_1)}$')
    if params_dict['press_scale'] == 0:
        ax_dict['deltap'].set( xlabel=r'$(P_2 - P_1) \cdot V \ [\varepsilon]$', ylabel=r'$q\Delta \varphi \ [\varepsilon]$')
    else:
        ax_dict['deltap'].set( xlabel=r'$(P_2 - P_1) \cdot V \ [10^{{{0[press_scale]}}} \cdot \varepsilon ]$'.format(params_dict), ylabel=r'$q\Delta \varphi \ [\varepsilon ]$')
    ax_dict['deltapT'].set( ylabel=r'$(P_2 - P_1) \cdot V \ [\varepsilon]$')


    # ax_dict['c'].set(xlim=[0, 1.2], xlabel='$T \ [\\varepsilon]$', ylabel='$\\Delta \\varphi^* \ [\\varepsilon]$')

    # K = 0.2
    # ax_dict['deltap'].plot(np.linspace(0,10,10), K*np.linspace(0,10,10), 'k--')
    # ax_dict['deltap'].plot(press_array, phi3_array_pure,'o',fillstyle='none', mec=c, markersize=13.,markeredgewidth=1.3,label=myfile)


    for ext in ['.png','.pdf']:
        for key, value in fig_dict.items():
            if key not in list_of_intern_axes:
            # fig_dict[key].show()
                fig_dict[key].savefig(key+ params_dict['out_suffix'] + ext)
    # fig_dict['deltap'].savefig('p_all.png')

if __name__ == '__main__':
    # main()


    def is_valid_file(parser, arg):
        """
        Check if arg is a valid file

        Args:
            parser (TYPE): Description
            arg (TYPE): Description

        Returns:
            TYPE: Description
        """
        arg = os.path.abspath(arg)
        if not os.path.exists(arg):
            parser.error("The file %s doesn't exist " % arg)
        else:
            return arg


    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('--lines', action='store_true', help='Plot the data with the lines style')

    parser.add_argument('-a', action='append', dest='myfiles_list',
                        type=str,
                        default=[],
                        help='Add values to analyze values to a list type (default: %(default)s)',
                        )
    parser.add_argument('--f', action='store_true', help='run test analysis')
    parser.add_argument('--c', action='store_true', help='run test analysis')
    parser.add_argument('--deltapT', action='store_true', help='run test analysis')
    parser.add_argument('--deltap', action='store_true', help='run test analysis')
    parser.add_argument('--Work', action='store_true', help='run test analysis')
    parser.add_argument('--press_normalize_ncount', action='store_true', help='run test analysis')

    parser.add_argument("-bs",
                        "--bound_suffix",
                        dest="bound_suffix",
                        default='interpolated',
                        type=str,
                        help="bound_suffix for the output csv file with the q1,q2 information (default: %(default)s)")


    parser.add_argument("-bp",
                        "--bound_prefix",
                        dest="bound_prefix",
                        default='Boundary',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")


    parser.add_argument("-cs",
                        "--cond_suffix",
                        dest="cond_suffix",
                        default='interpolated',
                        type=str,
                        help="cond_suffix for the output csv file with the q1,q2 information (default: %(default)s)")


    parser.add_argument("-cp",
                        "--cond_prefix",
                        dest="cond_prefix",
                        default='COND_q1_method',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")

    # add prefix and boundary
    # k_frames=9,startframe=2, endframe=30,trajskip=3

    parser.add_argument("-pl", "--plot_details", dest="plot_details",
                    default=1,
                    type=int,
                    help="how much detail do you want, 1 - just pure, 2 - pure Plus normalized by f1/f2, c1, c2; 3- normalized by f1/f2 c1/c2 and condensation, 4 - all together(default: %(default)s)")

    parser.add_argument("-ps", "--press_scale", dest="press_scale",
                    default=1,
                    type=int,
                    help="10^press_scale by which pressure will be scaled, 3 will mean that 5000 -> 5 [epsilon*10^3](default: %(default)s)")


    parser.add_argument("-os",
                        "--out_suffix",
                        dest="out_suffix",
                        default='_v1',
                        type=str,
                        help="suffix for the output plots so we don't overwrite them (default: %(default)s)")

     # parser.add_argument("-pn", "--press_normalize_ncount", dest="press_normalize_ncount",
     #                default=False,
     #                type=int,
     #                help="how much detail do you want, 1 - just pure, 2 - pure Plus normalized by f1/f2, c1, c2; 3- normalized by f1/f2 c1/c2 and condensation, 4 - all together(default: %(default)s)")
    parser.add_argument('--plot_pure', action='store_true', help='do you want to plot the pure raw data')


    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))


    main(args)
