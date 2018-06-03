#!/usr/bin/env python

# Program: Cluster3dFinal.py
# Purpose: Calculates clusters, using cluster labelling, analyzes their size
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python Cluster3dFinal.py -h for help,
# Requires: read_parameters.py -> to read trajectory, save_plots -> to vizualize the plot biggest cluster vs time
# AnalyzeCrystYamamoto -> to analyze crystallinity of the box

# THEORY:
# i don't have to vizualize everything
# i can just use the vizualizing by crystallinity within thershold
# in the User parameter of VMD
# and then perform the cluster analysis as a number analysis by using scipy.measurements
# the only thing is verification, how do I verify my results

# for every frame:
#     divide the box into a grid, by dividing
#     box -> planar section -> line section -> little box segment(dot)
#     then for every little box segment - we calculate crystallinity, and assign that crystallinity to w_cluster matrix
# after that we work with the w_cluster matrix -> we calculate biggest cluster in the frame,
# we calculate cluster size distribution, and plot a histogram

# PROBLEM: the user section is not updated dynamically
# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import scipy.ndimage
import MDAnalysis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import read_parameters
# import AnalyzeCrystYamamoto import get_bondlist_coords, get_eigvec, get_cryst
import AnalyzeCrystYamamoto
import save_plots
import get_path_names


def plot_hist(area, frame):
    """
    input: area -> sizes of the clusters, frame for animation
    plots histgoram
    this will plot the histogram
    and then convert it to gif
    PROBLEM : it can't create a histogram, beacuse all data is zero initialy
    """
    plt.ylim(0., 1.)
    plt.xlabel(r'$\mathrm{Cluster\ Size}$')
    plt.ylabel(r'$\mathrm{Population}$')
    plt.grid(True)
    plt.title(' Population =%f ' % (frame))
    n, bins, patches = plt.hist(area, 20,
                                range=(0, 100),
                                normed=1, facecolor='green',
                                alpha=0.75)
    np.savez('hist_%.5d' % int(frame), n, bins)
    plt.legend(loc='best')
    plt.savefig('hist_%.5d.png' % int(frame))
    plt.cla()


def get_biggest_cluster(w_cluster, frame):
    """
    input: w_cluster array
    output: the biggest cluster in the system
    method: uses the measurements module to analyze the clustering
    """
    # labeled_array, num_features = measurements.label(w_cluster)
    lw, num = scipy.ndimage.measurements.label(w_cluster)
    # lw += np.ones(len(lw),len(lw),len(lw))
    area = scipy.ndimage.measurements.sum(w_cluster, lw,
                                          index=np.arange(lw.max() + 1))
    areaImg = area[lw]
# the histogram of the data
    plot_hist(area, frame)
    return areaImg.max()


def create_grid_cryst(args):
    """
    input: args, obtained from read_traj_vmd
    output: timeseries of the crystallinity per frame
    method:
    for every frame:
        divide the box into a grid, by dividing
        box -> planar section -> line section -> little box segment(dot)
        then for every little box segment - we calculate crystallinity,
        and assign that crystallinity to w_cluster matrix
    after that we work with the w_cluster matrix
    -> we calculate biggest cluster in the frame,
    we calculate cluster size distribution, and plot a histogram
    """
    u = MDAnalysis.Universe(args.psffile, args.traj)
    a = u.selectAtoms("all")
    a.packIntoBox()
    psffile = get_path_names.get_filename(args.psffile)
    Nsub = args.nsub
# this is the gridding parameter,
# every dimenstion will be divided into Nsub cells
    cryst_all = 0.0
    box = u.trajectory.ts.dimensions[:-3]
    length_x = box[-1]
    grid_1d = np.linspace(0.0, length_x, Nsub+1, endpoint=True)
# grid in each dimension
    delta = grid_1d[1]-grid_1d[0]  # grid size

    timeseries = []  # this is for VMD userdata field
    sys_cryst = np.zeros(len(u.atoms))
# this is for VMD userdata field, it contains the userdata information
# on each atom per frame
# then we will put this info to timeseries :
# sys_cryst[ar_x.atoms.indices()] = box_cryst, timeseries.append(sys_cryst)

    time_frames = []  # this is for plot biggest_cluster_array vs time_frames
    bcluster_array = []
# this is for plot biggest_cluster_array vs time_frames
    w_cluster = np.zeros((Nsub, Nsub, Nsub))
# this is the w_cluster binary matrix, which will be used for cluster analysis

    # start dividing the simulation box into grid
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        # ar_z is a planar selection
        for k, z in enumerate(grid_1d[0:-1]):
            ar_z = u.selectAtoms("prop z >= {z0} and prop z < {z1}"
                                 .format(z0=str(z), z1=str(z+delta)))
            print ar_z.atoms
            print z
            # ater this step ar_y is a line
            for j, y in enumerate(grid_1d[0:-1]):
                ar_y = ar_z.selectAtoms("prop y >= {y0} and prop y < {y1}"
                                        .format(y0=str(y), y1=str(y+delta)))
                print ar_y.atoms
                # ater this step ar_x is a dot
                for i, x in enumerate(grid_1d[0:-1]):
                    ar_x = ar_y.selectAtoms("prop z >= {x0} and prop z < {x1}"
                                            .format(x0=str(x), x1=str(x+delta)))
                    print ar_x.atoms
                    bonds = AnalyzeCrystYamamoto.get_bondlist_coords(ar_x)
                    box_cryst = AnalyzeCrystYamamoto.get_cryst(bonds,
                         AnalyzeCrystYamamoto.get_eigvec(bonds,
                                                         threshold=args.threshold)

                    if box_cryst > args.threshold:
                        w_cluster[i, j, k] = 1.
                    elif 0. <= box_cryst <= args.threshold:
                        w_cluster[i, j, k] = 0.
                    else:
                        raise ValueError("box_cryst wrong type")
                    cryst_all += box_cryst
                    sys_cryst[ar_x.atoms.indices()] = box_cryst

        bcluster = get_biggest_cluster(w_cluster, ts.frame)
        # get the biggest cluster in the system
        # append information to arrays
        timeseries.append(sys_cryst)
        time_frames.append(ts.frame)
        bcluster_array.append(bcluster)
        cryst_all /= float(Nsub**3.0)
        print 'frame = {frame}, cryst = {cryst},\
              biggest_cluster = {bcluster}'\
              .format(frame=ts.frame,
                      cryst=cryst_all,
                      bcluster=bcluster)
        w_cluster = np.zeros((Nsub, Nsub, Nsub))
    # now we have arrays for all time frames,
    # lets either plot them, or put into vmd
    bcluster_array = np.array(bcluster_array)
    time_frames = np.array(time_frames)
    save_plots.save_plot(time_frames, bcluster_array,
                         name=psffile,
                         plotname="cluster",
                         xlabel="time",
                         ylabel="cluster size",
                         title="cluster"+psffile)
    return timeseries


def create_viz(args, timeseries):
    """
    input: args, from read_traj_vmd, timeseries: timeseries of the
    crystallinity per frame obtained from create_grid_cryst
    writes a vmd script, after initializing that it will create
    a user field data field
    we can put it into vmd and produce some vizualization
    """
    PSF = args.psffile
    DCD = args.traj
    # u = MDAnalysis.Universe(PSF, DCD)
    # psffile = get_path_names.get_filename(args.psffile)

    userdata = "box_cryst.txt"
    vmdscript = "viz_box.vmd"

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
                    $all set user $beta
                    puts "the data is loaded, frame = $frame"
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
    mol modselect 0 0 user > 0.3
    mol modselect 1 0 user < 0.25
    mol modcolor 1 0 ColorID 0
    # mol modstyle 1 0 CPK 1.000000 0.300000 12.000000 12.000000
    mol modstyle 1 0 Licorice 0.300000 12.000000 12.000000
    mol modmaterial 1 0 Diffuse

    mol modcolor 0 0 ColorID 0
    mol modcolor 0 0 ColorID 1
    # mol modstyle 0 0 Licorice 0.500000 46.000000 58.000000
    mol material 0 0 AOChalky
    mol modstyle 0 0 Licorice 0.600000 48.000000 58.000000
    mol modcolor 1 0 ColorID 7
    # mol modstyle 1 0 VDW 0.600000 27.000000

    # mol modselect 0 0 user > 0.6
    # mol modstyle 0 top Lines
    # mol modstyle 0 top VDW

    # mol modmaterial 0 0 AOChalky
    """.format(parameters)

    with open(vmdscript, 'w') as tcl:
        tcl.write(script+'\n')

    print("Wrote data trajectory {0} with res_cryst_array".format(userdata))
    print("Wrote VMD script {0}: 'source {0}' to load everything "
          .format(vmdscript))
    return None


def main():

    args = read_parameters.read_traj_vmd()
    timeseries = create_grid_cryst(args)
    # create_viz(args, timeseries)


if __name__ == '__main__':
    main()
