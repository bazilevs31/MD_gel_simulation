#!/usr/bin/env python

# Program: AnalyzeRMSD.py
# Purpose: calculates rmsd, by generating vmd file, which calculates it using vmd, then data is ploted
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeRMSD.py for help, 
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 


import numpy as np
from read_parameters import read_traj_vmd
from MDAnalysis import Universe
import os

#program for vizualizing crystallization

# PROBLEM: the user section is not updated dynamically - SOLVED

# PROBLEM: trajectory goes until the last frame, automatically, it doesn't know about args.endframe. - need to make a conditional 
# expression , so it analyzes the script 

def calc_rmsd(args):
    """
    input : args - args taken from read_traj_vmd
    it creates a Universe, 
    writes a vmdscript,
    executes it
    (vizualizing is done after by another function  plot_rdf )
    """

    psffile = os.path.splitext(args.psffile)[0]
    PSF = psffile+'.psf'
    DCD = args.traj
    # u = Universe(PSF, DCD)


    userdata = "gofr_"
    vmdscript = "calc_rdf.vmd"

    # write VMD loader script
    parameters = {'vmdfile': vmdscript,
                  'datfile': userdata,
                  'topology': PSF,
                  'trajectory': DCD,
                  'startframe': args.startframe,
                  'endframe': args.endframe,
                  'trajskip':args.trajskip}

    script = """\
    mol new "{0[topology]}"
    mol addfile "{0[trajectory]}" waitfor all

    # set outfile1 [open {0[datfile]} w]
    set sel [atomselect top "all"]
    set sel1 [atomselect top "all"]
    # set n [molinfo top get numframes]
    set outname "{0[datfile]}"
    set lastframe [molinfo top get numframes]

    for {{ set i 60}} {{ $i <  [ expr $lastframe - 3] }} {{ incr i 10}} {{

        $sel frame $i
        $sel update
        set gr0 [measure gofr $sel1 $sel delta 0.05 rmax 10.0 usepbc 1 selupdate 1 first $i last [expr $i + 2] step 1]
        set name [format "%s%d.dat" $outname $i]
        set outfile1 [open $name w]

        set r [lindex $gr0 0]
        set gr [lindex $gr0 1]
        foreach j $r k $gr {{
            puts $outfile1 [format "%.4f\t%.4f" $j $k]
        }}

    close $outfile1
    # $sel delete    
    }}
    exit
    """.format(parameters)


    with open(vmdscript, 'w') as tcl:
        tcl.write(script+'\n')

    # os.system("vmd -e {0[vmdfile]}".format(parameters))

    print("RDF data files will be written to  {0}_number.dat files".format(userdata))
    print("Wrote VMD script {0}: 'source {0}' to load everything ".format(vmdscript))
    return None

def plot_rmsd(args):
    """
    input : args - from read_traj_vmd
    to do this all the files in the directory will be listed, and saved to a list of files, 
    then each file in the list will be plotted and saved into a png file,
    then I will use the convert command of the imagemagick, to convert it to an animation .gif
    this function will plot the RDF
    """

    import os
    import glob
    import matplotlib.pyplot as plt

    psffile = os.path.splitext(args.psffile)[0]


    print 'Listing all profile/dat files'
    profilefilelist = glob.glob('gofr*.dat')
    print profilefilelist
    DATA=profilefilelist
    plt.ylim(0.,4.)
    print " data is red now reading the files"
    for i in DATA:
        # '''Read 2d array of z coordinates from file. Convert to float values
        # and wrap in a numpy array.'''
        with open(i) as f:
            data = [map(float, line.split()) for line in f]
        data = np.array(data)
        data = data.transpose()
        print "ploting file %r" % i
        plt.xlabel(r"distance $r$ in lj")
        plt.ylabel(r"radial distribution function $g(r)$")
        plt.title(str(i))
        plt.plot(data[0],data[1],linewidth=3)
        plt.savefig(i+'.png',figsize=(5,5),dpi=600)

    print " the files are plotted"
    print "creating gif animation"
    os.system("convert -delay 40 *.png animrdf.gif")
    return None

def main():
    """ the main function"""

    args = read_traj_vmd()
    psffile = os.path.splitext(args.psffile)[0]
    # calc_rdf(args)
    plot_rdf(psffile)

if __name__ == '__main__':
    main()