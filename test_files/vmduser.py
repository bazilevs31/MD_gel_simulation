#!/usr/bin/env python


# project a dynamic property on the structure in VMD using the User field

import numpy as np
from read_parameters import read_traj_vmd
import MDAnalysis
import MDAnalysis.analysis.align
from AnalyzeChain import get_bondlist_coords, get_chain_crystallinity
import os


args = read_traj_vmd()
psffile = os.path.splitext(args.psffile)[0]
PSF = psffile+'.psf'
DCD = args.traj
u = MDAnalysis.Universe(PSF, DCD)


CORE_selection = "all"
userdata = "adk_distance.dat"
vmdscript = "adk_distance.vmd"
timeseries = []

distances = np.ones(len(u.atoms))
# iterate through our trajectory
for ts in u.trajectory:
    # superimpose on the reference CORE (at t=0)
    # rmsd = MDAnalysis.analysis.align.alignto(u.atoms, ref.atoms, select=CORE_selection)
    # distances = np.sqrt(np.sum((u.atoms.positions - ref.atoms.positions)**2, axis=1))
    for res in u.residues:
                chords = get_bondlist_coords(res)
                res_g2 = get_chain_crystallinity(chords,args.neigh)
                distances[res.atoms.indices()] = res_g2


    timeseries.append(distances)

    print("Frame {0}: CORE RMSD before/after superposition: {1[0]:.1f} / {1[1]:.1f} A. "
          "min-max displacement: {2:.1f}...{3:.1f} A".format(ts.frame, rmsd, distances.min(), distances.max()))

# serialize: add a marker 'END' after each frame
marker = 'END'
with open(userdata, 'w') as data:
    for distances in timeseries:
        data.write("\n".join([str(x) for x in distances]))
        data.write("\n{0}\n".format(marker))

# write VMD loader script
parameters = {'datafile': userdata,
              'topology': PSF,
              'trajectory': DCD}

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
                $all set user $beta
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
mol addfile "{0[trajectory]}" waitfor all
loaduserdata "{0[datafile]}"
mol modcolor 0 top User
mol modstyle 0 top VDW
""".format(parameters)

with open(vmdscript, 'w') as tcl:
    tcl.write(script+'\n')

print("Wrote data trajectory {0} with distances".format(userdata))
print("Wrote VMD script {0}: 'source {0}' to load everything ".format(vmdscript))