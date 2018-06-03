#!/usr/bin/env python

import numpy as np
import re
import os
import argparse
import pickle

import pandas as pd
import Gel_dump_analysis

import matplotlib.pyplot as plt
plt.style.use('presentation')


def analyze_pressure_dump(temp_params, gel_params, press_params, folder, outname):
    """analyzes pressure and dumps the data"""
    fT = plt.figure(1)
    for gel_name, gel_temp_list, f1, f2 in zip(gel_params['Name'], gel_params['Temperatures'],gel_params['f1'],gel_params['f2']):
        print gel_name

        fT.clf()
        axT = fT.add_subplot(111)
        axT.set_xlim((0., 1.2))
        axT.set_xlabel(r'$T$')
        scale_power = 5
        scale_factor = 10.**scale_power
        axT.set_ylabel('$\Delta P (T) \cdot 10^{{-{0}}}$'.format(scale_power))
        axT.set_title(gel_name)
        # plt.grid('on')

        fz = plt.figure(2)

        phiarray = []
        temp_array = []
        density_phi_array = []
        legends = []
        press_dict={'p_ciKE':[],'p_gelELAS':[],'p_gelVIR':[],'p_gelPAIR':[]}
        error_press_dict={'p_ciKE':[],'p_gelELAS':[],'p_gelVIR':[],'p_gelPAIR':[]}

        for _t in gel_temp_list.split(','):
            _temp = float(_t)
            _color = temp_params.loc[temp_params['Temp'] == _temp]['Color'].values[0]
            _label = temp_params.loc[temp_params['Temp'] == _temp]['Label'].values[0]
            _tname = _t.replace('.', '')

            fz.clf()
            axz = fz.add_subplot(111)
            axz.set_title(_label)
            # axz.set_ylim((-0.7, 0.7))
            axz.set_xlabel('$z/l_z$')
            axz.set_ylabel('$P (z/l_z) \cdot 10^{{-{0}}}$'.format(scale_power))
            temp_array.append(_temp)
            print _temp, _color, _label
            for nm,m in zip(press_params['Pressure'],press_params['Color']):
                print "temp, tl, nm", gel_name, _temp, nm
                _labelp = press_params.loc[press_params['Pressure'] == nm]['Label'].values[0]
                df_means, df_stds = Gel_dump_analysis.analyze_pressure_dump("prof_"+gel_name+"_t"+_tname+"_ljeq1_coul_eq_coulsim_press")
                rel_press=scale_factor*df_means[nm]-scale_factor*df_means[nm][0]
                rel_error=scale_factor*2*df_stds[nm]
                axz.errorbar(df_means['Coord1'], rel_press,yerr=rel_error, fmt=m,label='${0}$'.format(_labelp))
                legends.append('${0}$'.format(_labelp))
                _fname = "./" + folder + "/abs"+gel_name+"_t"+_tname

                _p = np.average(rel_press[1:3]) - np.average(rel_press[-5:-1])
                _p = np.abs(_p)
                press_dict[nm].append(_p)
                error_press_dict[nm].append(rel_error.max())
            axz.legend(legends, loc='best')
            fz.savefig(_fname + '.png')
            pickle.dump(axz, file(_fname+'.pickle','w'))

        fT.tight_layout()
        legendsp = []
        for press_item, press_values in press_dict.iteritems():
            _colorp2 = press_params.loc[press_params['Pressure'] == press_item]['Color'].values[0]
            _markerp2 = press_params.loc[press_params['Pressure'] == press_item]['Marker'].values[0]
            _labelp2 = press_params.loc[press_params['Pressure'] == press_item]['Label'].values[0]
            # axT.plot(np.array(temp_array), np.array(press_values), _colorp2+_markerp2+'--',label=press_item)
            axT.errorbar(np.array(temp_array), np.array(press_values),yerr=error_press_dict[press_item], fmt=_colorp2+_markerp2+'--',label='${0}$'.format(_labelp2))
            legendsp.append('${0}$'.format(_labelp2))
            # np.savetxt(press_item + gel_name + '.csv', np.c_[temp_array, press_values])
            np.savetxt(press_item + gel_name+'test.csv', np.c_[temp_array,press_values], delimiter=',')


        temp_all_array = np.linspace(0,1.1,20)
        Ree = 50.
        deltaf = abs(float(f1) - float(f2))
        print deltaf
        press_all_array = 0.1*scale_factor * deltaf *100*temp_all_array/Ree**3
        # press_all_array =  0.1*0.4*100*temp_all_array/Ree**3
        print press_all_array
        # legendsp.append('${0}$'.format(_labelp2))
        # legendsp.append('IG theory')
        # box = axT.get_position()
        # axT.tight_layout
        # axT.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # axT.legend(legendsp, loc='center left', bbox_to_anchor=(1, 0.5))
        axT.legend(legendsp, loc='best', frameon=False)
        axT.plot(temp_all_array, press_all_array, 'r-')
        _fname2 = "./" + folder + "/" + outname + gel_name
        fT.savefig(_fname2 + '.png')
        pickle.dump(axT, file(_fname2+'.pickle','w'))

    return None

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

    analyze_pressure_dump(temp_params, gel_params, press_params, folder, outname)
