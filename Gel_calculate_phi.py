#!/usr/bin/env python

import numpy as np
import re
import os
import argparse
import pickle

import pandas as pd
import Gel_dump_analysis

# from Argparse import parser

import matplotlib.pyplot as plt
plt.style.use('presentation')

parser = argparse.ArgumentParser(description='input files in the csv format for analysis')

parser.add_argument('-t', '--temp', default='temperature_params.csv')
parser.add_argument('-p', '--press', default='press_params.csv')
parser.add_argument('-g', '--gel', default='gel_params.csv')
parser.add_argument('-o', '--out', default='press')
parser.add_argument('-fl', '--folder', default='pdf')
args = parser.parse_args()
folder = args.folder
outname = args.out

if not os.path.exists(folder):
    os.makedirs(folder)

temp_params = pd.read_csv(os.path.splitext(args.temp)[0]+'.csv')
gel_params = pd.read_csv(os.path.splitext(args.gel)[0]+'.csv')
press_params = pd.read_csv(os.path.splitext(args.press)[0]+'.csv')

"""analyzes pressure and dumps the data"""
fT = plt.figure(1)
legendsPhi = []

for gel_name, gel_temp_list in zip(gel_params['Name'], gel_params['Temperatures']):
    print gel_name

    fT.clf()
    axT = fT.add_subplot(111)
    axT.set_xlim((0., 1.2))
    axT.set_xlabel(r'$T$')
    axT.set_ylabel(r'$\Delta P (T) \cdot 10^{-7}$')
    axT.set_title(gel_name)
    # plt.grid('on')

    fz = plt.figure(2)

    phiarray = []
    temp_array = []
    density_nernst_array = []
    density_nernst_array_error = []
    legendsT = []

    for _t in gel_temp_list.split(','):
        _temp = float(_t)
        _color = temp_params.loc[temp_params['Temp'] == _temp]['Color'].values[0]
        _label = temp_params.loc[temp_params['Temp'] == _temp]['Label'].values[0]
        _tname = _t.replace('.', '')

        fz.clf()
        axz = fz.add_subplot(111)
        axz.set_title(_label)
        axz.set_ylim((-1, 1))
        axz.set_xlabel('$z$')
        axz.set_ylabel(r'$e \phi(z)/k_b$')


        f = gel_name + '_t' + _t
        data = np.loadtxt(f + '.dat')
        xx, yy = data[:,0], data[:,1]
        indices = np.where((xx>50.) & (xx<700))
        x = xx[indices]
        y = yy[indices]/500

        axz.plot(x, y,_color + 'o')
        legendsT.append(_label)
        temp_array.append(_temp)
        ymax = y.max()
        ymin = y.min()
        phi = np.abs(ymax - ymin)
        phiarray.append(np.abs(phi))
        _fname = "./" + folder + "/abs"+gel_name+"_t"+_tname
        axz.legend(legendsT, loc='best')
        fz.savefig(_fname + '.png')
        pickle.dump(axz, file(_fname+'.pickle','w'))


    _colorphi2 = gel_params.loc[gel_params['Name'] == gel_name]['Color'].values[0]
    _markerphi2 = gel_params.loc[gel_params['Name'] == gel_name]['Marker'].values[0]
    _labelphi2 = gel_params.loc[gel_params['Name'] == gel_name]['Label'].values[0]
    # axT.plot(np.array(temp_array), np.array(press_values), _colorp2+_markerp2+'--',label=press_item)
    axT.plot(np.array(temp_array), np.array(phiarray), _colorphi2+_markerphi2+'--')
    legendsPhi.append('${0}$'.format(_labelphi2))
    axT.legend(legendsPhi, loc='best', fancybox=False, framealpha=1.)
    _fname2 = "./" + folder + "/" + outname + gel_name
    fT.savefig(_fname2 + '.png')
    pickle.dump(axT, file(_fname2+'.pickle','w'))


