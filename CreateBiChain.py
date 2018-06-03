#!/usr/bin/env python

# Program: CreateBiChain.py
# Purpose: create ndisperse melt
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreateBiChain.py -h for help,
# example: python CreateBiChain.py -bi 1 100 -bi 2 200 , 1*C100 + 2*C200
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version
# import os
import numpy as np
import read_parameters
import save_plots


def main():
    args = read_parameters.read_create_ndisp()
    bichains = np.array(args.bichains)
    n = bichains[:, 1]
    bins = bichains[:, 0]
    PDI = save_plots.save_create_hist(bins, n, type='bi')
    save_plots.save_write_to_file(n, bins, PDI)

if __name__ == '__main__':
    main()

# modify it
# maybe use this :
# numpy.random.choice
# create a range of chain lengths , then select them randomly with certain
# weights
# for example
# create a np.range(40,160):every 10 then select random stuff ,with probability that is another exponentially decreasing function that will represent a poisson distribution

