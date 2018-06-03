#!/usr/bin/env python

# Program: CreatePDChain.py
# Purpose: create polydisperse chains systems
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreatePDChain.py for help,
# example: python CreatePDChain.py -mu 100 -s 20 -nch 100
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import read_parameters
import save_plots
import matplotlib.pyplot as plt

# modify it
# maybe use this :
# numpy.random.choice
# create a range of chain lengths , then select them randomly with certain
# weights
# for example
# create a np.range(40,160):every 10 then select random stuff ,with probability that is another exponentially decreasing function that will represent a poisson distribution


def create_hist(args):
    """
    histogramms data and return bin values and population of bins
    input: args , from read_parameters
    output: bins, n
    laplace distribution
    x = np.arange(-8., 8., .01)
>>> pdf = np.exp(-abs(x-loc)/scale)/(2.*scale)
    """
    mu = args.mu
    sigma = args.mu
    Nchains = args.Nchains
    Nbins = args.Nbins

    x = np.random.laplace(mu, sigma, Nchains)
    # x_hist = np.linspace(20, 200, Nbins)
    x_hist = np.arange(5, 1000, Nbins)
    n, bins = np.histogram(x, bins=x_hist)
    # n, bins, patches = plt.hist(x, bins=x_hist, facecolor='green', alpha=0.75)
    # plt.show()
    return bins[:-1], n


def main():
    args = read_parameters.read_create_pd()
    bins, n = create_hist(args)
    PDI = save_plots.save_create_hist(bins, n,
                                      type='pd',
                                      mu=args.mu,
                                      sigma=args.sigma,
                                      Nchains=args.Nchains)
    save_plots.save_write_to_file(n, bins, PDI)

if __name__ == '__main__':
    main()

# for fitting data
    # from pylab import *
    # from numpy import loadtxt
    # from scipy.optimize import leastsq

    # fitfunc  = lambda p, x: p[0]*exp(-0.5*((x-p[1])/p[2])**2)+p[3]
    # errfunc  = lambda p, x, y: (y - fitfunc(p, x))

    # filename = "gaussdata.csv"
    # data     = loadtxt(filename,skiprows=1,delimiter=',')
    # xdata    = data[:,0]
    # ydata    = data[:,1]

    # init  = [1.0, 0.5, 0.5, 0.5]

    # out   = leastsq( errfunc, init, args=(xdata, ydata))
    # c = out[0]

    # print "A exp[-0.5((x-mu)/sigma)^2] + k "
    # print "Parent Coefficients:"
    # print "1.000, 0.200, 0.300, 0.625"
    # print "Fit Coefficients:"
    # print c[0],c[1],abs(c[2]),c[3]

    # plot(xdata, fitfunc(c, xdata))
    # plot(xdata, ydata)

    # title(r'$A = %.3f\  \mu = %.3f\  \sigma = %.3f\ k = %.3f $' %(c[0],c[1],abs(c[2]),c[3]));

    # show()
