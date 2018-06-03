#!/usr/bin/env python
import MDAnalysis
import numpy as np
import cellgrid
import numba
import math
import scipy
import scipy.ndimage
import test_cell as analyze_cell


def test_1():
    """
    runs a test to figure out how well the cluster algorithms work
    In [12]: areaImg
    Out[12]:
    array([[[ 2.,  2.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]],

        [[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]],

        [[ 3.,  3.,  3.],
            [ 0.,  0.,  0.],
            [ 0.,  1.,  0.]]])
    """

    cryst_box = np.zeros((3,3,3),dtype=np.int)
    cryst_box[0,0,:-1] = 1
    cryst_box[2,0,:] = 1
    cryst_box[2,2,1]= 1


    lw, num = scipy.ndimage.measurements.label(cryst_box)
        # lw += np.ones(len(lw),len(lw),len(lw))
    area = scipy.ndimage.measurements.sum(cryst_box, lw,
                                            index=np.arange(lw.max() + 1))

    areaImg = area[lw]
    x, y, z= areaImg.nonzero()
    analyze_cell.plot_cells_3d(x, y, z, areaImg[x,y,z], 1, name='test3d')
    num_clusters_array = len(((area>0).nonzero())[0])
    num_large_clusters_array = (len(((area>10).nonzero())[0]))
    analyze_cell.plot_hist(area, 1, name='testplot')


def test_2():
    """
    runs a test to figure out how well the cluster algorithms work
    """
    N = 10
    cryst_box =  np.random.choice(np.arange(2), size=(N,N,N), p = [0.8, 0.2])
    lw, num = scipy.ndimage.measurements.label(cryst_box)
        # lw += np.ones(len(lw),len(lw),len(lw))
    area = scipy.ndimage.measurements.sum(cryst_box, lw,
                                            index=np.arange(lw.max() + 1))

    areaImg = area[lw]
    x, y, z= areaImg.nonzero()
    analyze_cell.plot_cells_3d(x, y, z, areaImg[x,y,z], 1, name='test3d')
    num_clusters_array = len(((area>0).nonzero())[0])
    num_large_clusters_array = (len(((area>10).nonzero())[0]))
    analyze_cell.plot_hist(area, 1, name='testplot')
    print num_clusters_array
    print num_large_clusters_array

