#!/usr/bin/env python
import MDAnalysis
import numpy as np
import cellgrid
import numba
import math
import scipy
import scipy.ndimage
import read_parameters
import AnalyzeLog
import matplotlib
matplotlib.use('Agg')

FIG_X, FIG_Y = 8, 8
FIG_DPI = 80



@numba.jit(nopython=True)
def calc_whether_aligned(loc_bonds, loc_align_vec, threshold=0.85,
                         cryst_perc_threshold=0.5):
    """
    calculates number of aligned parts in subdomain for Yamamoto Crystallinity of Data, XDEF2
    :param bonds: np.array(N,3, dtype=np.float32)
    :param align_vec: np.array(3, dtype=np.float32)
    :param threshold:  float, default .99
    :return: int number aligned parts

    .. warning:: both $bonds$ and $align_vec$ should be normalized
    """
    cos = 0.
    k = 0
    tot_bonds = loc_bonds.shape[0]
    for i in range(tot_bonds):
        cos = 0.
        for j in range(3):
            cos += loc_align_vec[j] * loc_bonds[i, j]
        if math.fabs(cos) > threshold:
            k += 1
    if (k/tot_bonds>cryst_perc_threshold):
        return 1
    else:
        return 0

def vmd_clusters(args, timeseries):
    """
    puts the values of the cluster size into the vmd readable format
    then writes the .vmd file
    """
    marker = 'END'
    userdata = "tmp_viz.txt"
    vmdscript = "viz_dat.vmd"
    PSF = 'mdc40.psf'
    DCD = 'traj.dcd'
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



# @numba.jit(nopython=True)
def get_eigvec(bonds, loc_qab, loc_align_vec):
    """
    calculates the eigen vector that show the alignment directions of chains
    :param bonds: np.array((N,3), dtype=np.float32) of normalized
    bonds vectors of a selected domain
    :param loc_qab: np.array((3,3), dtype=np.float32) just for effectivity,
    you actually don't need it, but this way faster
    :param loc_align_vec: np.array(3,dtype=np.float32) eigen vector,
     which corresponds to orientation of
    selected bonds
    :return: loc_align_vec
    """
    nq = bonds.shape[0]
    loc_qab = 1.5*(1./nq)*np.einsum('ij,ik->jk', bonds, bonds) - 0.5*np.eye(3,dtype=np.float32)
    vals, vecs = np.linalg.eig(loc_qab)
    loc_align_vec = vecs[:, 0]
    return loc_align_vec

def plot_cells_3d(x_local, y_local, z_local, color_array, frame, name):
    """docstring for plot_cells_3d"""
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(FIG_X, FIG_Y), dpi=FIG_DPI)
    ax = fig.add_subplot(111, projection='3d')
    # marker_size = FIG_X*FIG_DPI/8.
    marker_size = 150
    ax.set_title('time ' + str(frame))
    ax.scatter(x_local, y_local, z_local, marker='s',c=color_array, s=marker_size, zdir='z',
               cmap=plt.cm.rainbow)
    plt.savefig(name+'_%f.png' % frame)


def calc_clusters(args):
    """
    calculates clusters
    """
    psf_name = args.psffile
    dcd_name = args.traj
    # psf_name = './test_mono_100/mono100n1000.psf'
    # dcd_name = './test_mono_100/skipquenchsim.dcd'
    u = MDAnalysis.Universe(psf_name, dcd_name)
    box = u.trajectory.ts.dimensions
    timeseries = []

    # getting the dumpskip parameters
    if args.filedumpskip:
        dumpskip = AnalyzeLog.get_dumpskip(filename='dumpskip.txt',
                                           keyparameter="dump" +
                                           args.keyparameter+"skip")
        # dumpskip = args.dumpskip
        timestep = AnalyzeLog.get_timestep(filename='dumpskip.txt',
                                           keyparameter="time" +
                                           args.keyparameter
                                           + "step")
    else:
        print "using given parameters"
        dumpskip = args.dumpskip
        timestep = args.timestep

    # decide the skipping information
    if args.auto_traj_skip:
        trajskip = int(u.trajectory.numframes/120.)
        if trajskip < 1:
            trajskip = 1
    elif args.auto_traj_skip is False:
        trajskip = args.trajskip
    else:
        raise ValueError("trajskip needs to be provided")


    qab = np.empty((3, 3), dtype=np.float32)
    align_vec = np.empty(3, dtype=np.float32)

    cell_size = 3.0
    frames_array = []
    num_clusters_array = []
    num_large_clusters_array = []

    angle_threshold=0.8
    perc_threshold=0.3
    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        # print "doing loop", ts.frame
        # current time in 10^6 LJ tau
        cur_time = ts.frame*dumpskip*timestep/1000000.
        u.SYSTEM.pack_into_box()
        coords = (u.bonds.atom2.positions + u.bonds.atom1.positions)*0.5
        bonds = u.bonds.atom2.positions - u.bonds.atom1.positions
        norm = np.linalg.norm(bonds, axis=1)
        bonds /= norm[:, None]
        cg = cellgrid.CellGrid(box[:3], float(cell_size), coords)
        cryst_box = np.zeros(cg._ncells, dtype=np.int)
        cryst_index_array = np.zeros(cg._total_cells)
        record_array = np.zeros(len(u.atoms), dtype=np.float32)

        # print cg


        for cell in cg:
                cell_indices = cell.indices
                align_vec = get_eigvec(bonds[cell_indices], qab, align_vec)
                k = calc_whether_aligned(bonds[cell_indices],
                                        align_vec,
                                        threshold=angle_threshold,
                                        cryst_perc_threshold=perc_threshold)
                cryst_index_array[cell.index] = k
                cryst_box[cell.address] = k

        lw, num = scipy.ndimage.measurements.label(cryst_box)
            # lw += np.ones(len(lw),len(lw),len(lw))
        area = scipy.ndimage.measurements.sum(cryst_box, lw,
                                                index=np.arange(lw.max() + 1))
        # print "lw", lw
        areaImg = area[lw]
        x, y, z = areaImg.nonzero()
        # print cell
        # print cell.index
        for cell in cg:
            record_array[cell.indices] = areaImg[cell.address]
        # print record_array
        # print record_array.shape
        # print record_array.nonzero()
        timeseries.append(record_array)
        # print cell.address
        # print cryst_box.shape
        # print areaImg.shape
        # areaImg - connects the information about the cluster
        # size in terms of number of cells cluster entake

        # now I need to put all the atoms of that cell this index

        # serialize: add a marker 'END' after each frame
    userdata = "tmp_viz.txt"
    vmdscript = "viz_dat.vmd"
    marker = 'END'
    with open(userdata, 'w') as data:
        for res_cryst_array in timeseries:
            data.write("\n".join([str(x) for x in res_cryst_array]))
            data.write("\n{0}\n".format(marker))

    # write VMD loader script
    parameters = {'vmdfile': vmdscript,
                  'datafile': userdata,
                  'topology': psf_name,
                  'trajectory': dcd_name,
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

        # timeseries.append(res_cryst_array)

        # now I need to use areaImg
        # to program the thing


        #initially all atoms are zeros
        # here I assign atoms the color of their cluster size
        # here I assign


    #     # print "frame = ", ts.frame, " number clusters of 3 ", len((areaImg>2).nonzero())
    #     print "time mln LJ tau = ", cur_time, " number clusters of 3 ", len(((area>0).nonzero())[0])
    #     plot_cells_3d(x, y, z, areaImg[x,y,z], cur_time, name='demo')
    #     frames_array.append(cur_time)
    #     num_clusters_array.append(len(((area>0).nonzero())[0]))
    #     num_large_clusters_array.append(len(((area>10).nonzero())[0]))
    #     print area
    #     plot_hist(area, cur_time, name='hist')

    #     plot_graph(np.array(frames_array),
    #             np.array(num_clusters_array),
    #             np.array(num_large_clusters_array),
    #             ts.frame,
    #             threshold=angle_threshold,
    #             cryst_perc_threshold=perc_threshold,
    #             name='plot')
    np.savetxt(u.filename+'.csv',np.c_[np.array(frames_array),
                np.array(num_clusters_array),
                np.array(num_large_clusters_array)], delimiter=',')

def main():
    args = read_parameters.read_traj_vmd()
    calc_clusters(args)

if __name__ == '__main__':
    main()
