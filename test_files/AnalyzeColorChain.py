#!/usr/bin/env python

# Program: AnalyzeColorChain.py
# Purpose: calculates, crystallinity, writes a vmd file, which will vizualize the crystallinity
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python AnalyzeColorChain.py for help, 
# Requires: read_parameters.py, AnalyzeChain.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version 

import numpy as np
from read_parameters import read_traj_vmd
from MDAnalysis import Universe
from AnalyzeChain import get_bondlist_coords, get_chain_crystallinity
import os

#program for vizualizing crystallization

# PROBLEM: the user section is not updated dynamically
def main():

    args = read_traj_vmd()
    psffile = os.path.splitext(args.psffile)[0]
    PSF = psffile+'.psf'
    DCD = args.traj
    u = Universe(PSF, DCD)


    userdata = "tmp_viz.txt"
    vmdscript = "viz_dat.vmd"
    timeseries = []

    res_cryst_array = np.zeros(len(u.atoms))
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        for res in u.residues:
                    chords = get_bondlist_coords(res)
                    res_g2 = get_chain_crystallinity(chords,args.threshold,args.neigh)    #calculate crystallinity of a chain
                    # print res.atoms.indices()
                    res_cryst_array[res.atoms.indices()] = res_g2
        timeseries.append(res_cryst_array)

        print "frame %d" % ts.frame
        print "result ", np.average(res_cryst_array)

    # serialize: add a marker 'END' after each frame
    marker = 'END'
    with open(userdata, 'w') as data:
        for res_cryst_array in timeseries:
            data.write("\n".join([str(x) for x in res_cryst_array]))
            data.write("\n{0}\n".format(marker))

    # write VMD loader script
    parameters = {'vmdfile': vmdscript,
                  'datafile': userdata,
                  'topology': PSF,
                  'trajectory': DCD,
                  'startframe': args.startframe,
                  'endframe': args.endframe,
                  'trajskip': int(args.trajskip)}

    script = """\
    proc loaduserdata { fname } {
        set all [atomselect top all]
        set frame 0
        set data [open $fname r]
        while { [gets $data line] != -1 } {
            set value [string trim $line]
            switch -- [string range $value 0 2] {
                END {
                    $all frame $frame
                    # $all set user $beta
                    $all set vx $beta
                    set beta {}
                    incr frame
                }
                default {
                    lappend beta $line
                }
            }
        }
    }
    """ + """
    mol new "{0[topology]}"
    animate read dcd "{0[trajectory]}" skip {0[trajskip]}  waitfor all
    # mol addfile "{0[trajectory]}" type {{dcd}} first 0 last -2 step 1 waitfor
    # mol addfile "{0[trajectory]}" waitfor all 
    animate goto 0
    loaduserdata "{0[datafile]}"
    color change rgb  0 0.1 0.2 0.7 ;# blue
    color change rgb  1 0.7 0.2 0.1 ;# red
    color change rgb  3 0.7 0.4 0.0 ;# orange
    color change rgb  4 0.8 0.7 0.1 ;# yellow
    color change rgb  7 0.1 0.7 0.2 ;# green
    color change rgb 10 0.1 0.7 0.8 ;# cyan
    color change rgb 11 0.6 0.1 0.6 ;# purple
    # color Display Background white #black
    color Display Background white
    # mol selupdate 0 0 1
    # mol colupdate 0 0 1
    mol modcolor 0 top User
    mol material AOChalky
    mol addrep 0
    mol modselect 0 0 user > 0.27
    mol modselect 1 0 user < 0.20
    mol modcolor 1 0 ColorID 0
    # mol modstyle 1 0 CPK 1.000000 0.300000 12.000000 12.000000
    mol modstyle 1 0 Licorice 0.300000 12.000000 12.000000
    mol modmaterial 1 0 Diffuse

    mol modcolor 0 0 ColorID 0
    mol modcolor 0 0 ColorID 1
    # mol modstyle 0 0 Licorice 0.500000 46.000000 58.000000
    # mol material 0 0 AOChalky
    mol modmaterial 0 0 AOChalky
    mol modstyle 0 0 Licorice 0.600000 48.000000 58.000000
    mol modcolor 1 0 ColorID 7
    # mol modstyle 1 0 VDW 0.600000 27.000000

    # mol modselect 0 0 user > 0.6
    # mol modstyle 0 top Lines
    # mol modstyle 0 top VDW

    # mol modmaterial 0 0 AOChalky

    box_molecule top


    """.format(parameters)

    with open(vmdscript, 'w') as tcl:
        tcl.write(script+'\n')

    print("Wrote data trajectory {0} with res_cryst_array".format(userdata))
    print("Wrote VMD script {0}: 'source {0}' to load everything ".format(vmdscript))


if __name__ == '__main__':
    main()