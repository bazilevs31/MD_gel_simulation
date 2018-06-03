#!/usr/bin/env python

# Program: AnalyzeGelInit.py
# Purpose: analyzes initial gel structure
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeGelInit.py -h for help,

# Copyright (c) 2016 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='program to create nodes, counter ions and salt ions (both cations and anions)')
parser.add_argument('-f', '--filename', default='mol')
parser.add_argument('-n', '--nbins', default=2)
args = parser.parse_args()
input_data = (args.filename).rsplit('.', 1)[0]

print args

FILE = input_data
PSF = FILE + ".psf"
DCD = FILE + ".dcd"
DATA = FILE + ".data"


try:
    u = MDAnalysis.Universe(PSF, DCD)
except Exception as e:
    u = MDAnalysis.Universe(DATA)
# u = MDAnalysis.Universe(PSF, DCD)

u.universe.atoms.pack_into_box()
counter_ions = u.select_atoms("name 3")
nodes = u.select_atoms("name is 1")
gel = u.select_atoms("name is 1 or name is 2")
x_countions = counter_ions.positions[:,0]
x_gel = gel.positions[:,0]
Nbins = int(args.nbins)

hist_array, bin_all = np.histogram(x_countions, Nbins)
print bin_all
print hist_array
plt.plot(bin_all[:-1], hist_array, 'go-', label='counter ions',markersize=20)
hist_array, bin_all = np.histogram(x_gel, Nbins)
print bin_all
print hist_array

print " =======       ========="
print " ======= total ========="
print " =======       ========="

Nci = counter_ions.n_atoms
Nnodes = nodes.n_atoms
Ngel = gel.n_atoms
Nmons =  Ngel - Nnodes
fraction = Nci/float(Ngel)
print "counter ions", Nci
print "nodes", Nnodes
print "gel", Ngel
print "fraction", Nci/float(Ngel)

plt.title('gel={gel}, Cions={cions}, fraction={f}'.format(gel=Ngel,cions=Nci, f=float(fraction)))

plt.plot(bin_all[:-1], hist_array, 'ro-', label='gel',markersize=10)
plt.legend()

plt.savefig(FILE+'.png')
