#!/usr/bin/env python

# Program: AnalyzeBondACF.py
# Purpose: calculate bond autocorrelation function
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeBondACF.py -h for help,
# example: python AnalyzeBondACF.py
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version
# import os
import numpy as np
import MDAnalysis
import read_parameters
import AnalyzeChain
import save_plots
import get_path_names
import numba


# @numba.jit
def calc_ci(bonds, C):
    Nbonds = C.shape[0]
    count = 0.*C
    for beg_bond in range(0, Nbonds):
        for ui in range(beg_bond, Nbonds):
            C[ui] += np.dot(bonds[beg_bond], bonds[ui])
            count[ui] += 1
    return C/count


def calc_bondacf(args):
    """
    calculates the stem length by calculating bond correlation functino
    bondacf: bond correlation function analysis for L_p, when this
                parameter is provided, then program analyzes first points
                and fits them, to understand persistance length.
                C(n) = <u_0 u_n >
    method: create an array C,
    loop trajectory:
        loop residues:
            loop bonds:
                analyze angle between first bond vector and i
    average results
    plot using save_plots
    """
    psffile = get_path_names.get_filename(args.psffile)
    u = MDAnalysis.Universe(psffile + '.psf', args.traj)
    Nbonds = len(u.residues[0])-1
    Nres = len(u.residues)
    print "Nres = %r" % Nres
    print "Nbonds = %r" % Nbonds
    bonds = np.array((Nbonds, 3))
    C = np.zeros(Nbonds)
    k = 0
    # for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
    # print "time = %r" % ts.frame
    for res in u.residues:
        bonds = AnalyzeChain.get_bondlist_coords(res)
        k += 1
        C = calc_ci(bonds, C)
    # C /= float(k)
    C /= len(u.residues)
    print C
    ui_array = np.arange(Nbonds)
    save_plots.save_plot(ui_array, C,
                         plotname='stemlength',
                         xlabel='n',
                         ylabel='bond')


def main():
    args = read_parameters.read_traj_vmd()
    calc_bondacf(args)

if __name__ == '__main__':
    main()
