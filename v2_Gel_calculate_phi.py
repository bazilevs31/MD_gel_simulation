#!/usr/bin/env python

import numpy as np
import re
import os
import argparse
import pickle

import pandas as pd
import Gel_dump_analysis
from gridData import Grid

import pint
ureg = pint.UnitRegistry()

from lmfit.models import StepModel, LinearModel

# from argparse import parser

import matplotlib.pyplot as plt
plt.style.use('presentation')


def analyze_dx(f):
    """analyzes dx file"""
    g = Grid(f + ".dx")
    nx, ny, nz = g.grid.shape
    Lx = g.edges[0][0]
    Ly = g.edges[1][0]
    Lzl = g.edges[2][0]
    Lzr = g.edges[2][-1]
    # g_all = g.grid[nx//2][ny//2]

    # g_vmd = g_all
    # L = 100.

    g_all = np.zeros(nz)

    count = 0

    for i in range(nx):
        for j in range(ny):
            g_all += g.grid[i][j]
            count += 1
    r_array = np.linspace(Lzl, Lzr, nz)
    # vmd_unit = 25*ureg.mV
    # vmd_unit = 1.
    # units = nu.e/(25*nu.mV)
    # print units
    yy = g_all[2:-2]/count
    return r_array, yy

def fit_data(xdata, ydata, xfull):
    """takes the data
    returns yfull, amplitude"""

    step_mod = StepModel(form='logistic', prefix='step_')
    line_mod = LinearModel(prefix='line_')

    pars =  line_mod.make_params(intercept=ydata.min(), slope=0)
    pars += step_mod.guess(ydata, x=xdata, center=250, amplitude=1.)

    # pars.add('amp_1', value=0.5, min=0, max=1)
    # pars.add('step_amplitude',value=1.0, min = 0.1, max=1.2)
    mod = step_mod + line_mod
    out = mod.fit(ydata, pars, x=xdata)

    print(out.fit_report())

    # plt.plot(x_array, out.init_fit, 'k--')
    # plt.plot(x, out.best_fit, 'r-')

    # step_mod.func

    yfull = step_mod.func(xfull, amplitude=out.best_values['step_amplitude'], center=out.best_values['step_center'], sigma=out.best_values['step_sigma'], form=out.best_values['step_form'])
# step_mod.func(x, amplitude=1.0, center=0.0, sigma=1.0, form='linear')
    yfull2 = line_mod.func(xfull, slope=out.best_values['line_slope'], intercept=out.best_values['line_intercept'])

    return out.best_values['step_amplitude'], yfull+yfull2

def get_delta_stepfit(_x, _y, a=100., b=700.):
    """calculates phi as a difference between roots of the fitting derivative"""
    indices = np.where((_x>a) & (_x<b))
    x = _x[indices]
    y = _y[indices]
    xfull = np.linspace(x.min(), x.max(), 1000)
    phi, yfull = fit_data(x, y, xfull)
    return np.abs(phi)

def get_delta_polyfit(_x, _y, a=100., b=700.):
    """calculates phi as a difference between roots of the fitting derivative"""
    indices = np.where((_x>a) & (_x<b))
    x = _x[indices]
    y = _y[indices]

    z = np.polyfit(x, y, 4)
    p = np.poly1d(z)
    p2 = np.polyder(p)
    x_roots = np.roots(p2)
    y_roots = p(x_roots)
    delta_phi = np.abs(y_roots.max() - y_roots.min())
    return delta_phi

def get_delta_minmax(_x, _y, a=100., b=700.):
    """calculates phi as a difference between the max and min value"""
    indices = np.where((_x>a) & (_x<b))
    x = _x[indices]
    y = _y[indices]
    ymax = y.max()
    ymin = y.min()
    _phi = np.abs(ymax - ymin)
    return np.abs(_phi)

def analyze_phi(temp_params, gel_params, press_params, folder, outname, method='minmax'):
    """
    analyzes potential

    """
    fT = plt.figure(1)
    fT.clf()
    axT = fT.add_subplot(111)
    # axT.set_xlim((0., 1.2))
    axT.set_ylim((0., 1.2))
    axT.set_xlabel(r'$T$')
    axT.set_ylabel(r'$e\Delta \varphi (T)/k_b$')
    legendsPhi = []

    fz = plt.figure(2)
    axz = fz.add_subplot(111)
    axz.set_ylim((-1, 1))
    axz.set_xlabel('$z$')
    axz.set_ylabel(r'$e \phi(z)/k_b$')

    for gel_name, gel_temp_list in zip(gel_params['Name'], gel_params['Temperatures']):
        print gel_name

        # plt.grid('on')


        phiarray = []
        temp_array = []
        temp_all_array = np.linspace(0.,1.1, 20)
        density_nernst_array = []
        density_nernst_array_error = []
        legendsT = []


        for _t in gel_temp_list.split(','):
            print "now doing T= ", _t
            _temp = float(_t)
            _color = temp_params.loc[temp_params['Temp'] == _temp]['Color'].values[0]
            _label = temp_params.loc[temp_params['Temp'] == _temp]['Label'].values[0]
            _tname = _t.replace('.', '')
            f = gel_name + '_t' + _tname
            axz.set_title(f)
            # data = np.loadtxt(f + '.dat')
            # xx, yy = data[:,0], data[:,1]
            xx, yy = analyze_dx(f)
            xx, yy = xx, yy/500.
            # print xx,yy

            axz.plot(xx, yy, _color + 'o--')

            legendsT.append(_label)
            temp_array.append(_temp)
            if method == 'polyfit':
                phi=get_delta_polyfit(xx,yy)
                print "method", method, phi
            elif method == 'minmax':
                phi=get_delta_minmax(xx,yy)
                print "method", method, phi
            elif method == 'stepfit':
                phi=get_delta_stepfit(xx,yy)
                print "method", method, phi
            phiarray.append(phi)
            _fname = "./" + folder + "/abs"+gel_name+"_t"+_tname
        axz.legend(legendsT, loc='best')
        fz.savefig(_fname + '.png')
        # plt.gcf().clear()
        axz.cla()
        pickle.dump(axz, file(_fname+'.pickle','w'))
        # axz.clear()


        _colorphi2 = gel_params.loc[gel_params['Name'] == gel_name]['Color'].values[0]
        _labelphi2 = gel_params.loc[gel_params['Name'] == gel_name]['Label'].values[0]
        _f1 = gel_params.loc[gel_params['Name'] == gel_name]['f1'].values[0]
        _f2 = gel_params.loc[gel_params['Name'] == gel_name]['f2'].values[0]
        # axT.plot(np.array(temp_array), np.array(press_values), _colorp2+_markerp2+'--',label=press_item)
        print gel_name, phiarray
        axT.plot(np.array(temp_array), np.array(phiarray), _colorphi2+'o')
        legendsPhi.append('MD ${0}$'.format(_labelphi2))
        axT.plot(temp_all_array, -np.log(float(_f1)/float(_f2))*temp_all_array, _colorphi2+'--')
        legendsPhi.append('Nernst ${0}$'.format(_labelphi2))
    axT.legend(legendsPhi, loc='best', fancybox=False, framealpha=1.)
    _fname2 = "./" + folder + "/" + outname + gel_name
    fT.savefig(_fname2 + '.png')
    pickle.dump(axT, file(_fname2+'.pickle','w'))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='input files in the csv format for analysis')
    parser.add_argument('-t', '--temp', default='temperature_params.csv')
    parser.add_argument('-p', '--press', default='press_params.csv')
    parser.add_argument('-g', '--gel', default='gel_params.csv')
    parser.add_argument('-o', '--out', default='press')
    parser.add_argument('-fl', '--folder', default='pdf')
    args = parser.parse_args()
    folder = args.folder

    if not os.path.exists(folder):
        os.makedirs(folder)

    outname = args.out

    temp_params = pd.read_csv(os.path.splitext(args.temp)[0]+'.csv')
    gel_params = pd.read_csv(os.path.splitext(args.gel)[0]+'.csv')
    press_params = pd.read_csv(os.path.splitext(args.press)[0]+'.csv')

    _method='minmax'
    analyze_phi(temp_params, gel_params, press_params, folder, _method+outname, method=_method)

    # for _method in ['polyfit', 'minmax','stepfit']:
    #     analyze_phi(temp_params, gel_params, press_params, folder, _method+outname, method=_method)
