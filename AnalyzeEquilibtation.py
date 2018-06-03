#!/usr/bin/env python

# Program: AnalyzeEquilibration.py
# Purpose: using all files listed below analyze equilibration and plot the data
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeEquilibration.py -h for help, 
# Requires: read_parameters.py, AnalyzeEnd2End, AnalyzeLog, AnalyzeRDF, AnalyzeRgRe

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 

from read_parameters import read_traj_vmd, read_log
import os
from AnalyzeEnd2End import calc_r2n
from AnalyzeLog import analyze_log
from AnalyzeRgRe import calc_rgre,calc_rmsd
from AnalyzeRDF import calc_rdf,plot_rdf

# in this file i will Use 

# AnalyzeEnd2End.py for calculating R(n)/n and ploting

# AnalyzeRDF.py for calculating RDF, -> using VMD, gnuplot

# AnalyzeRMSD.py for calculating RMSD - > not implemented yet , use VMD for calculating rmsd, save to file,

# AnalyzeRgRe.py for calcuting Rg(t), Re(t)

# AnalyzeLog.py for ploting the thermodynamics information


def main():
    """
    analyze equilibration of melt 
    """

    args = read_traj_vmd()
    print " "
    print " ====================== "
    print " done inputing the trajectory"
    print args
    print ""
    print " ====================== "
    print " "


    calc_r2n(args)
    print " "
    print " ====================== "
    print " done r2(n)/n"
    print ""
    print " ====================== "
    print " "
    calc_rgre(args)
    print " "
    print " ====================== "
    print " done rgre"
    print ""
    print " ====================== "
    print " "

    calc_rmsd(args)
    print " "
    print " ====================== "
    print " done rmsd"
    print ""
    print " ====================== "
    print " "
    calc_rdf(args)
    plot_rdf(args)
    print " "
    print " ====================== "
    print " done rdf"
    print ""
    print " ====================== "
    print " "
    # logargs = read_log()
    # print " "
    # print " ====================== "
    # print " done inputing the log"
    # print logargs
    # print ""
    # print " ====================== "
    # print " "
    # analyze_log(logargs)
    return None

if __name__ == '__main__':
    main()