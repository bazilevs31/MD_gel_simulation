#!/usr/bin/env python

# Program: AnalyzeLog.py
# Purpose: calculates all thermodynamic data, in a single file in a format Value_i(Step)
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeLog.py -h for help,
# Requires: read_parameters.py, save_plot_log, pizlog

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import os
from log import log as pizlog
from termcolor import colored
# import log as pizlog
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse
import termcolor

def read_log():
    """
    provide information for logfile
    output: args , containing logfile, name of the output file
    """

    parser = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-l", "--log", dest="logfile",
                        default="./figures/log.kremer",
                        help="read logfile to analyze")

    parser.add_argument("-s",
                    "--stride",
                    dest="initoffset",
                    default=1,
                    type=int,
                    help="How many first steps are to be skipped, because they are ususally too big")
    parser.add_argument("-nv",
                    "--Nevery",
                    dest="Nevery",
                    default=10,
                    type=int,
                    help="Plot every this step")

    parser.add_argument("-o",
                    "--output",
                    dest="outfile",
                    default="sim",
                    help="results of the logfile thermo info ploting will be saved here")

    parser.add_argument("-al",
                    "--all",
                    action="store_true",
                    dest="all",
                    default=False,
                    help="---Plot all---: do you want to plot all log files in ./figures(default: %(default)s)")
    ######### log plotting

    parser.add_argument("-lp",
                    "--logplot",
                    action="store_true",
                    dest="logplot",
                    default=False,
                    help="Do you want time(or N) to be plotted with log(t) scale? (default: %(default)s)")

    args = parser.parse_args()

    return args

def save_plot_log(x, thermodata, **kwargs):
    """
    plots data from thermo log file
    """

    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    parameters = dict()
    input_dict = {key: value for (key, value) in kwargs.items()}
    Nevery = int(input_dict['Nevery'])
    parameters['plotname'] = input_dict.get('plotname', 'log')
    parameters['name'] = input_dict.get('name', 'poly')
    parameters['xlabel'] = input_dict.get('xlabel', 'time, 10^6 lj units')
    parameters['figuresdir'] = figuresdir
    print termcolor.colored('parameters have been red', 'green')
    print parameters
    with matplotlib.backends.backend_pdf.PdfPages('{0[plotname]}{0[name]}.pdf'.
                                                  format(parameters),
                                                  'w') as pdf:

        for i, element in enumerate(thermodata):
            print element
            plt.figure(figsize=(10, 6))
            plt.plot(x[::Nevery], thermodata[str(element)][::Nevery], 'b-')
            plt.title(str(element))
            plt.ylabel(str(element))
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
        d = pdf.infodict()
        d['Title'] = 'equilibration of polymer thermodynamic log'
        d['Author'] = u'Vasiliy M. Triandafilidi'
        d['Subject'] = 'How to plot lammps log file thermodynamic entries'
        d['Keywords'] = 'equilibration lammps log pizza'

    plt.figure()
    time_array = x[::Nevery]
    Temp_array = thermodata['Temp'][::Nevery]
    Vol_array = thermodata['Volume'][::Nevery]
    PotEng_array = thermodata['PotEng'][::Nevery]
    plt.ylim(0.3, 1.1)
    plt.title('Thermodynamic Parameters Evoution')
    plt.ylabel('Normalized Thermodynamic parameters')
    plt.xlabel('${0[xlabel]}$'.format(parameters))
    plt.plot(time_array, Temp_array, 'r-',
     label=r'$Temp$')
    plt.plot(time_array, Vol_array/Vol_array.max(), 'g-',
     label=r'$\frac{vol}{vol_{max}}$')
    plt.plot(time_array, PotEng_array/PotEng_array.max(), 'b-',
     label=r'$\frac{E_{pot}}{E_{potmax}}$')
    # plt.legend(loc='best')
    # plt.legend()
    plt.legend(loc='best')
    plt.savefig(r'publication{0[name]}.pdf'.format(parameters))
    np.savez('{0[figuresdir]}/{0[name]}'.format(parameters), time_array, Temp_array, Vol_array)
    return None


def plot_func(x, y, **kwargs):
    """
    plots stuff
    input: x, y, **kwargs
    x,y : arrays of the same length to plot
    kwargs : additional parameters to make the plot better
        plotname : name of algorithm used for calculating x,y
        (for example ICC,YCryst)
            ICC : Individual Chain Crystallinity Parameter
            YC : Yamamoto Crystallinity Parameter
            Cos : Average cos2*8eta
            Rg : array of Rg crystallinity parameters
            Re : array of Re : end to end distance evolution
            Bead mean:square displacement
            Sq : static structure factor
            Lml : lamellae analysis
            bondacf: bond correlation function analysis for L_p, when this
                parameter is provided, then program analyzes first points
                and fits them, to understand persistance length.
                C(n) = <u_0 u_n >

        name : name of the psffile
        xlog, ylog : whether you want log scale or not
        title : name of the plot in the title
        xlabel, ylabel : labels
        duplicate: True, False
            whether you want to duplicate results of the simulations in
            ~/results_all/ folder

    TODO: bondacf to work properly

    """
    x = np.asarray(x, dtype='float32')
    y = np.asarray(y, dtype='float32')
    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    parameters = dict()
    input_dict = {key: value for (key, value) in kwargs.items()}
    parameters['plotname'] = input_dict.get('plotname', 'deformation')
    parameters['name'] = input_dict.get('name', 'poly')
    parameters['xlabel'] = input_dict.get('xlabel', '$\varepsilon')
    parameters['ylabel'] = input_dict.get('ylabel', 'pressure')
    parameters['title'] = input_dict.get('title', 'deformation')
    parameters['xlog'] = input_dict.get('xlog', False)
    parameters['ylog'] = input_dict.get('ylog', False)
    parameters['ymin'] = input_dict.get('ymin', 0.0)
    parameters['ymax'] = input_dict.get('ymax', 2.0)
    parameters['duplicate'] = input_dict.get('duplicate', False)
    parameters['figuresdir'] = figuresdir
    print termcolor.colored('parameters have been red', 'green')
    print parameters

    plt.ylabel('${0[ylabel]}$'.format(parameters))
    plt.xlabel('${0[xlabel]}$'.format(parameters))
    plt.title('{0[title]}'.format(parameters))
    plt.grid(True)

    plt.plot(x, y, 'bo-', lw=1.5)
    np.savez('{0[plotname]}{0[name]}'.format(parameters), x, y)
    plt.savefig('{0[plotname]}{0[name]}.pdf'
                .format(parameters))

def analyze_log(args):
    """
    input : args -> from read_log
    this function performs analysis of MD trajectory thermodata log

    # this program analyzes equilibration log
    # to analyze it we use pizza.py program
    # parameters to analyze
    # Etot(t) T(t)
    # Epot(t) V(t)
    # Ebond(t) msd(t)
    #
    # after analyzing it
    # plot all quantities using subplot
    #
    """
    namelog = args.logfile
    initoffset = args.initoffset
    outfile = args.outfile
    Nevery = args.Nevery
    if args.logplot is True:
        print "xlogscale"
        print "{log(time),\  lj\ units}}"
    else:
        print "regular xscale"
        print "{(time),\  lj\ units}}"

    # if we have filelogskip flag set to True then we get the logskip information from a file,
    # else we use the one that comes from -ls flag

        #get timestep as well, and keyparameter

    ###### ---------- getting the names for files-----------------

    head, outfile= os.path.split(outfile)  # stripping the name so only name of the file is left
    outfile = os.path.splitext(outfile)[0]

    head, tail = os.path.split(namelog)
    name = os.path.splitext(tail)[0]

    ###### ---------main part--------------------

    lg = pizlog(namelog) #pizza.py part we read the thermodynamic info into lg variable
    thermodata = dict()
    #now we read this lg to create a python dictionary that is easy to use  and to plot
    for names in lg.names:
        print "reading %s " % names
        thermodata[str(names)]=np.array(lg.get(str(names)))
        thermodata[str(names)]=np.delete(thermodata[str(names)], range(initoffset))
    logtime=thermodata['Time'] # getting time, and normalizing it to the 10^6 lj units
    # logtime=thermodata['Step']*logskip*timestep/factor # getting time, and normalizing it to the 10^6 lj units

    # print logtime

    ###### ---------plotting  part--------------------
    print " all data is set for plotting"
    # save_plot_log(logtime,thermodata, outfile, name, Nevery, logplot=args.logplot)
    # plot_func(thermodata['Lx'], thermodata['Press'])
    save_plot_log(logtime, thermodata,
                             plotname='log',
                             name=name,
                             logplot=args.logplot,
                             Nevery=Nevery,
                             xlabel='time, lj\ units')

    return None
def main():
    args = read_log()
    analyze_log(args)

if __name__ == '__main__':
    main()
