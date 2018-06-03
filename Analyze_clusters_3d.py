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
    np.savetxt(name+'_t_' + str(frame)+'.csv', np.c_[x_local, y_local,
     z_local,color_array], delimiter=',')
    plt.savefig(name+'_%f.png' % frame)


def plot_hist(cluster_array, frame, name='hist'):
    """docstring for plot_cells_3d"""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(FIG_X, FIG_Y), dpi=FIG_DPI)
    ax = fig.add_subplot(111)
    ax.set_title('time = ' + str(frame))
    ax.set_ylim(0,40)
    # ax.set_xlim(0,20)
    # N, bins, patches = ax.hist(cluster_array, bins=5, normed=1)
    if np.all(np.array(cluster_array)==0.):
        pass
    else:
        print "cluster array", cluster_array
        N, bins, patches = ax.hist(cluster_array, bins=5)
        N_patches = len(patches)
        jet = plt.get_cmap('rainbow', N_patches)
        for i in range(N_patches):
            patches[i].set_facecolor(jet(i))
        # plt.colorbar.ColorbarBase(ax, cmap=jet,
                                # orientation='horizontal')
    plt.savefig(name+'%f.png' % frame)

def plot_graph(x, y_all, y_large, frame,
               threshold,
               cryst_perc_threshold,
               name='plot'):
    """
    plots number of crystallites vs time
    for parameters
    cryst_perc_threshold
    threshold
    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(FIG_X, FIG_Y), dpi=FIG_DPI)
    ax = fig.add_subplot(111)
    ax.set_ylim(0,40)
    ax.set_xlim(0,0.6)
    ax.set_title("{0} time {1}_params:_th_{2}perc_th{3}".format(name,
                                                                         frame,
                                                                         threshold,
                                                                         cryst_perc_threshold))
    ax.plot(x, y_all, 'go-',label='all')
    ax.plot(x, y_large, 'ro-',label='>30')
    ax.legend(loc='best')
    plt.savefig(name+str(frame)+'.png')


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
        x, y, z= areaImg.nonzero()
        # print "frame = ", ts.frame, " number clusters of 3 ", len((areaImg>2).nonzero())
        print "time mln LJ tau = ", cur_time, " number clusters of 3 ", len(((area>0).nonzero())[0])
        plot_cells_3d(x, y, z, areaImg[x,y,z], cur_time, name='demo')
        frames_array.append(cur_time)
        num_clusters_array.append(len(((area>0).nonzero())[0]))
        num_large_clusters_array.append(len(((area>10).nonzero())[0]))
        print area
        plot_hist(area, cur_time, name='hist')

        plot_graph(np.array(frames_array),
                np.array(num_clusters_array),
                np.array(num_large_clusters_array),
                ts.frame,
                threshold=angle_threshold,
                cryst_perc_threshold=perc_threshold,
                name='plot')
    np.savetxt(u.filename+'.csv',np.c_[np.array(frames_array),
                np.array(num_clusters_array),
                np.array(num_large_clusters_array)], delimiter=',')

def main():
    args = read_parameters.read_traj_vmd()
    calc_clusters(args)

if __name__ == '__main__':
    main()
