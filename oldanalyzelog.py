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
from read_parameters import read_log
from log import log as pizlog
from termcolor import colored
import re
import save_plots
# import log as pizlog
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

def get_timestep(filename='dumpskip.txt',keyparameter='logquenchskip'):
    """searches the keyparameter in the text file and returns dumpskip if it exists"""
    print colored('keyparameter = %s' % keyparameter, 'cyan')
    lines = filter(None, (line.rstrip() for line in open(filename)))
    flag = False
    print "looking for dumpskip keyword = %s" % keyparameter
    for line in lines:
        print "analyzing file %s"  % line
        if flag==False:
            if re.search('(\S+) (\d+.\d+)', line).group(1)==keyparameter:
                dumpskip = re.search('(\S+) (\d+.\d+)', line).group(2)
                flag = True
    if flag==False:
        raise ValueError("the keyparameter = %s is not found" % keyparameter)
    print colored('parameters have been red', 'cyan')
    print colored('timestep = %f' % float(dumpskip), 'green')
    return float(dumpskip)

def get_dumpskip(filename='dumpskip.txt',keyparameter='logquenchskip'):
    """searches the keyparameter in the text file and returns dumpskip if it exists"""
    print colored('keyparameter = %s' % keyparameter, 'cyan')
    lines = filter(None, (line.rstrip() for line in open(filename)))
    flag = False
    print "looking for dumpskip keyword = %s" % keyparameter
    for line in lines:
        print "analyzing file %s"  % line
        if flag==False:
            if re.search('(\S+) (\d+)', line).group(1)==keyparameter:
                dumpskip = re.search('(\S+) (\d+)', line).group(2)
                flag = True
    if flag==False:
        raise ValueError("the keyparameter = %s is not found" % keyparameter)
    print colored('parameters have been red', 'cyan')
    print colored('dumpskip = %d' % int(dumpskip), 'green')
    return int(dumpskip)


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

    initoffset = args.initoffset
    outfile = args.outfile
    Nevery = args.Nevery
    factor = 1e6 # the time units
    if args.logplot is True:
        print "xlogscale"
        print "{log(time),\ 10^6 lj\ units}}"
    else:
        print "regular xscale"
        print "{(time),\ 10^6 lj\ units}}"

    # if we have filelogskip flag set to True then we get the logskip information from a file,
    # else we use the one that comes from -ls flag
    if args.filelogskip is True:
        logskip = get_dumpskip(filename='dumpskip.txt',keyparameter="log"+args.keyparameter+"skip")
        # dumpskip = get_dumpskip(filename='dumpskip.txt',keyparameter=args.keyparameter)
        timestep = get_timestep(filename='dumpskip.txt',keyparameter="time"+args.keyparameter+"step")

        #get timestep as well, and keyparameter
    else:
        print "using logskip and timestep from argparse"
        logskip = args.logskip
        timestep = args.timestep

    namelog = args.logfile
    print " I am doing %s " % namelog

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
    print timestep
    logtime=thermodata['Step']*timestep/factor # getting time, and normalizing it to the 10^6 lj units
    # logtime=thermodata['Step']*logskip*timestep/factor # getting time, and normalizing it to the 10^6 lj units

    # print logtime

    ###### ---------plotting  part--------------------
    print " all data is set for plotting"
    # save_plot_log(logtime,thermodata, outfile, name, Nevery, logplot=args.logplot)
    save_plots.save_plot_log(logtime, thermodata,
                             plotname='log',
                             name=name,
                             logplot=args.logplot,
                             dumpskip=logskip,
                             timestep=timestep,
                             Nevery=Nevery,
                             xlabel='time,\ 10^6 lj\ units')

    return None
def main():
    args = read_log()
    analyze_log(args)

if __name__ == '__main__':
    main()
