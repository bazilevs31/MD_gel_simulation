#!/usr/bin/env python


import os
import glob
import sys
import subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from numpy import *
from numpy.random import rand
from pylab import pcolor, show, colorbar, xticks, yticks
from pylab import *


# def PlotRDF():
#     """

#     this function takes arguments of which files to plot and plots them into .png files, after that 
#     it converts all of the png files into a single .gif 

#     """

print 'Listing all profile/dat files'
profilefilelist = glob.glob('gofr*.dat')
raw_input('Press ENTER to continue...')
print profilefilelist
DATA=profilefilelist


for i in DATA:
    def get_data(fname=i):
        '''Read 2d array of z coordinates from file. Convert to float values
        and wrap in a numpy array.'''
        with open(fname) as f:
            data = [map(float, line.split()) for line in f]
            return np.array(data)
for i in DATA:    
    def my_plot(x, y):
        fig = plt.figure()
        matplotlib.rc('font', size=14)
        matplotlib.rc('figure', figsize=(5, 4))
        pylab.clf()
        pylab.plot(radii, rdf, linewidth=3)
        pylab.xlabel(r"distance $r$ in $\AA$")
        pylab.ylabel(r"radial distribution function $g(r)$")
        pylab.savefig(outfile + ".pdf")
        pylab.savefig(outfile + ".png")
        print "Figure written to ./figures/" + filename  + "{pdf,png}"

        ax.plot(x,y)
        # ax = fig.gca(projection='3d')
        # ax.plot_surface(x, y, z, rstride=5, cstride=5,cmap="binary",linewidth=0.1)
        # ax.set_zlim3d(0.0,4.0)  
        ax.set_xlabel('X',fontsize=16,fontweight="bold")
        ax.set_ylabel('Y',fontsize=16,fontweight="bold")
        # ax.set_zlabel('h(X,T)',fontsize=16,fontweight="bold")
        # plt.show()
       savefig(os.getcwd()+DATA+'.pdf',figsize=(5,5),dpi=600)
#        savefig(os.getcwd()+DATA+'.pdf',figsize=(5,5),dpi=600)    



if __name__ == '__main__':
            z = get_data()
            z = z.transpose()
            my_plot(z[0],z[1])