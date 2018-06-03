#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.interpolate
import matplotlib.pyplot as plt
import os
import termcolor
import argparse

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def interopolate_df_boundary(df1, plot_flag=False, **kwargs):
    """interpolates and extrapolates data"""

    # newdf = pd.DataFrame(columns=['T','q1','q2'])
    xdata = df1['T']
    y1 = df1['z1']
    y2 = df1['z2']
    y3 = df1['dl']


    nans, x= nan_helper(y1)
    y1[nans] = np.interp(x(nans), x(~nans), y1[~nans])
    y2[nans] = np.interp(x(nans), x(~nans), y2[~nans])
    y3[nans] = np.interp(x(nans), x(~nans), y3[~nans])

    f1 = scipy.interpolate.interp1d(xdata, y1,fill_value="extrapolate")
    f2 = scipy.interpolate.interp1d(xdata, y2,fill_value="extrapolate")
    f3 = scipy.interpolate.interp1d(xdata, y3,fill_value="extrapolate")

    temprange = np.linspace(0.1,1.1,11)
    y1new = f1(temprange)   # use interpolation function returned by `scipy.interpolate.interp1d`
    y2new = f2(temprange)   # use interpolation function returned by `scipy.interpolate.interp1d`
    y3new = f3(temprange)   # use interpolation function returned by `scipy.interpolate.interp1d`
    # newdf = pd.DataFrame([temprange, y1new, y2new],columns=['T','q1','q2'])
    _df = pd.DataFrame([temprange, y1new, y2new, y3new])
    df2 = _df.T
    df2.columns = ['T','z1','z2','dl']

    if plot_flag:
        plt.clf()
        plt.title('{0[myfile]}'.format(kwargs))
        plt.plot(df1['T'], df1['dl'], 'go')
        # plt.plot(df2['T'], df2['z2'], 'r.--')
        plt.plot(df2['T'], df2['dl'], 'r.--')
        plt.xlabel('T')
        plt.ylabel('$\Delta l$')
        plt.savefig('{0[csv_in_str]}_{0[myfile]}_{0[csv_out_suffix]}.png'.format(kwargs))
        # plt.show()
    return df2

def interopolate_df(df1, plot_flag=False, **kwargs):
    """scipy.interpolates and extrapolates data"""

    # newdf = pd.DataFrame(columns=['T','q1','q2'])
    xdata = df1['T']
    y1 = df1['q1']
    y2 = df1['q2']
    # y3 = df1['param']
    print(y1)

    nans, x= nan_helper(y1)
    y1[nans] = np.interp(x(nans), x(~nans), y1[~nans])
    y2[nans] = np.interp(x(nans), x(~nans), y2[~nans])
    # y3[nans] = np.interp(x(nans), x(~nans), y3[~nans])

    f1 = scipy.interpolate.interp1d(xdata, y1,fill_value="extrapolate")
    f2 = scipy.interpolate.interp1d(xdata, y2,fill_value="extrapolate")
    # f3 = scipy.interpolate.interp1d(xdata, y3,fill_value="extrapolate")

    temprange = np.linspace(0.1,1.1,11)
    y1new = f1(temprange)   # use interpolation function returned by `interp1d`
    y2new = f2(temprange)   # use interpolation function returned by `interp1d`
    # y3new = f3(temprange)   # use interpolation function returned by `interp1d`
    # newdf = pd.DataFrame([temprange, y1new, y2new],columns=['T','q1','q2'])
    _df = pd.DataFrame([temprange, y1new, y2new])
    df2 = _df.T
    df2.columns = ['T','q1','q2']
    df2['param'] = (1-df2['q1'])/(1-df2['q2'])

    if plot_flag:
        plt.plot(df1['T'], df1['q1'], 'go')
        plt.plot(df2['T'], df2['q1'], 'r.--')
        plt.savefig(df2.to_csv('{0[csv_in_str]}_{0[myfile]}_{0[csv_out_suffix]}.csv'.format(kwargs)))
        # plt.show()
    return df2

def run_analysis(filelist, csv_in_str='COND_v1', csv_out_suffix='interpolated'):
    """loops through the list of files, interpolates the data,
    if needed and saves a modified csv file"""

    for myfile in filelist:
        csvname = '{csv_in_str}_{myfile}.csv'.format(csv_in_str=csv_in_str, myfile=myfile)
        if os.path.exists(csvname):
            print('file exists')
        # df = pd.read_csv('COND_q_'+myfile+'.csv')
            df = pd.read_csv(csvname)
            mydf = df.drop(columns=['Unnamed: 0'])
            if 'Bound' in csv_in_str:
                print("Boundary is in csv_in_str")
                newdf = interopolate_df_boundary(mydf, plot_flag=True, csv_in_str='COND_v1', myfile=myfile, csv_out_suffix='interpolated')
            elif ('COND' in csv_in_str) or ('Cond' in csv_in_str):
                print("Boundary is NOT in csv_in_str")
                print(csvname)
                newdf = interopolate_df(mydf, plot_flag=True, csv_in_str='COND_v1', myfile=myfile, csv_out_suffix='interpolated')
            newdf.to_csv('{csv_in_str}_{myfile}_{csv_out_suffix}.csv'.format(csv_in_str=csv_in_str, myfile=myfile, csv_out_suffix=csv_out_suffix))
        else:
            print('FILE doesnt exist')
    return None

    # it doesn't do the interpolation properly but we can use it for now


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a', action='append', dest='files',
                        type=str,
                        default=[],
                        help='Add values to analyze values to a list type (default: %(default)s)')

    parser.add_argument("-cs",
                        "--csv_prefix",
                        dest="csv_prefix",
                        default='COND_v1',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")
    parser.add_argument("-name",
                        "--out_suffix",
                        dest="out_suffix",
                        default='interpolated',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")
    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))

    # _tmp1 = 'sc1844_nm100_l50_fl01_fr1'
    # _tmp2 = 'sc1844_nm100_l50_fl01_fr02'
    # _tmp3 = 'sc1844_nm100_l50_fl025_fr05'
    # _tmp4 = 'sc1844_nm100_l50_fl01_fr05'
    # filelist = [_tmp1,_tmp2,_tmp3,_tmp4]

    run_analysis(args.files, csv_in_str=args.csv_prefix, csv_out_suffix=args.out_suffix)
