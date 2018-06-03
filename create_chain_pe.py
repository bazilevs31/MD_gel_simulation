#!/usr/bin/env python
"""Summary
"""

# Program: create_chain_v2.py
# Purpose: creates a chain between given nodes
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python create_chain_v2.py -h for help,

# Copyright (c) 2016 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import math
# given two points A and B
# create random points that are adjacent to A and B
# interpolate between them repeat until the distance is 1
import random
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import distances

def dist(x1, y1, x2, y2):
    """Summary
    calculates distance between two points
    Args:
        x1 (float): point
        y1 (float): point
        x2 (float): point
        y2 (float): point

    Returns:
        float: distance
    """
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

def dist_3d(x1, y1,z1, x2, y2,z2):
    """Summary

    Args:
        x1 (float): point
        y1 (float): point
        z1 (float): point
        x2 (float): point
        y2 (float): point
        z2 (float): point

    Returns:
        float: dist
    """
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


def build_chain_3d_final(x1,x2, N, b, bounds,verbose=False):
    """
    build a 3d chain between given points using self avoiding walk
    and interpolation
    Args:
        x1 (np array x3): array of the atom coordinates
        x2 (np array x3): array of the atom coordinates
        N (int): polymer length
        b (float): bond length
        bounds (np array x3): box
        verbose (bool, optional): print or

    Returns:
        (np array (N-2,3)): array of the mons coordinates
        (float): bond length
    """
    if verbose:
        print "size is N", N
    x2 = distances.get_pbc(x1, x2, bounds)
    # N += 2
    # N += 1
    xx1, yy1, zz1, xx2, yy2, zz2 = x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]
    # if N*b <= dist_3d(xx1, yy1, zz1, xx2, yy2, zz2):
    #     print "error should be longer"
    # else:
    X = np.zeros(N)
    Y = np.zeros(N)
    Z = np.zeros(N)


    # ax.scatter(xx1, yy1, zz1, 'ro')
    # ax.scatter(xx2, yy2, zz2, 'ro')
    for i in range(N/2):
        if verbose:
            print "doing iter i ", i
        X[i] = xx1
        Y[i] = yy1
        Z[i] = zz1
        X[-(i+1)] = xx2
        Y[-(i+1)] = yy2
        Z[-(i+1)] = zz2

        alphai_1 = math.pi*random.random()
        alphai_2 = math.pi*random.random()
        # theta=sin-1(z/R), alpha = phi
        # x=Rcos(theta)cos(phi), y=Rcos(theta)sin(phi), z=Rsin(theta)
        rand_z_1 = random.random()
        rand_z_2 = random.random()
        theta_1 = math.asin(rand_z_1/b)
        theta_2 = math.asin(rand_z_2/b)

        xx1_i = xx1 + b * math.cos(alphai_1) * math.cos(theta_1)
        yy1_i = yy1 + b * math.sin(alphai_1) * math.cos(theta_1)
        zz1_i = zz1 + b * math.sin(theta_1)

        xx2_i = xx2 + b * math.cos(alphai_2) * math.cos(theta_2)
        yy2_i = yy2 + b * math.sin(alphai_2) * math.cos(theta_2)
        zz2_i = zz2 + b * math.sin(theta_2)

        d = dist_3d(xx1_i, yy1_i, zz1_i, xx2_i, yy2_i, zz2_i)
        Nleft = N - 2 * (i+1)
        dr_left = Nleft * b
        if verbose:
            print "nleft ", Nleft
        if d < dr_left:
            if verbose:
                print "normal step, dist", d, " dist left = ", dr_left, " 1:x,y ", xx1_i, yy1_i, " 2:x,y ", xx2_i, yy2_i
            xx1, yy1, zz1 = xx1_i, yy1_i, zz1_i
            xx2, yy2, zz2 = xx2_i, yy2_i, zz2_i
        elif d >= dr_left:
            # do some magic with reiterating the last step
            if verbose:
                print "too close line step, dist", d, " dist left = ", dr_left
            # xx, yy1 = xx_i, yy1_i
            # x2, y2 = x2_i, y2_i
            x = np.linspace(xx1, xx2, Nleft, endpoint=False, retstep=False)
            y = np.linspace(yy1, yy2, Nleft, endpoint=False, retstep=False)
            z = np.linspace(zz1, zz2, Nleft, endpoint=False, retstep=False)
            for j in range(Nleft):
                X[i+j+1], Y[i+j+1], Z[i+j+1] = x[j], y[j], z[j]
                # X[i+j], Y[i+j], Z[i+j] = x[j], y[j], z[j]
                # print "steps ", x[j], y[j], z[j]
            break
        else:
            print "error"


    # print X
    # print Y
    # print Z
    # return (X, Y,Z)
    if verbose:
        print "size result is ", len(X[1:-1])
    return np.vstack((X[1:-1], Y[1:-1],Z[1:-1])).T, b
    # return np.vstack((X[1:], Y[1:],Z[1:])).T, b


def build_chain_3d(x1, y1,z1, x2, y2,z2, b, N):
    """Summary

    Args:
        x1 (TYPE): Description
        y1 (TYPE): Description
        z1 (TYPE): Description
        x2 (TYPE): Description
        y2 (TYPE): Description
        z2 (TYPE): Description
        b (TYPE): Description
        N (TYPE): Description

    Returns:
        TYPE: Description
    """
    print "size is N", N
    verbose = False
    if N*b <= dist_3d(x1, y1, x2, y2):
        print "error should be longer"
    else:
        X = np.empty(N)
        Y = np.empty(N)
        Z = np.linspace(z1, z2, N, endpoint=False, retstep=False)


        for i in range(N/2):
            X[i] = x1
            Y[i] = y1
            X[-(i+1)] = x2
            Y[-(i+1)] = y2

            alphai_1 = math.pi*random.random()
            alphai_2 = math.pi*random.random()

            x1_i = x1 + b * math.cos(alphai_1)
            y1_i = y1 + b * math.sin(alphai_1)

            x2_i = x2 + b * math.cos(alphai_2)
            y2_i = y2 + b * math.sin(alphai_2)

            d = dist_3d(x1_i, y1_i, x2_i, y2_i)
            Nleft = N - 2 * (i+1)
            dr_left = Nleft * b
            if d < dr_left:
                if verbose:
                    print "normal step, dist", d, " dist left = ", dr_left, " 1:x,y ", x1_i, y1_i, " 2:x,y ", x2_i, y2_i
                x1, y1 = x1_i, y1_i
                x2, y2 = x2_i, y2_i
            elif np.allclose(d, dr_left):
                # connect with a straight line
                if verbose:
                    print "straight line step, dist", d, " dist left = ", dr_left
                x1, y1 = x1_i, y1_i
                x2, y2 = x2_i, y2_i
                for j in range(Nleft):
                    x1_i = x1 + b * i
                    y1_i = y1 + b * i
                    X[i+j+1], Y[i+j+1] = x1_i, y1_i
                    if verbose:
                        print "steps ", x1_i, y1_i
                break
            elif d > dr_left:
                # do some magic with reiterating the last step
                if verbose:
                    print "too close line step, dist", d, " dist left = ", dr_left
                # x1, y1 = x1_i, y1_i
                # x2, y2 = x2_i, y2_i
                x = np.linspace(x1, x2, Nleft, endpoint=False, retstep=False)
                y = np.linspace(y1, y2, Nleft, endpoint=False, retstep=False)
                for j in range(Nleft):
                    X[i+j+1], Y[i+j+1] = x[j], y[j]
                    if verbose:
                        print "steps ", x[j], y[j]
                break
            else:
                print "error"


        print X
        print Y
        print Z
        # return np.vstack((X, Y,Z)), b
        return np.hstack((X, Y,Z)), b
        # plt.plot(X, Y, 'bo-')

def build_chain(x1, y1, x2, y2, b, N):
    """Summary

    Args:
        x1 (TYPE): Description
        y1 (TYPE): Description
        x2 (TYPE): Description
        y2 (TYPE): Description
        b (TYPE): Description
        N (TYPE): Description

    Returns:
        TYPE: Description
    """
    plt.xlim(-10,10)
    plt.ylim(-10,10)

    if N*b <= dist(x1, y1, x2, y2):
        print "error should be longer"
    else:
        X = np.empty(N)
        Y = np.empty(N)


        plt.plot(x1, y1, 'ro',markersize=20)
        plt.plot(x2, y2, 'ro',markersize=20)
        for i in range(N/2):
            X[i] = x1
            Y[i] = y1
            X[-(i+1)] = x2
            Y[-(i+1)] = y2

            alphai_1 = math.pi*random.random()
            alphai_2 = math.pi*random.random()

            x1_i = x1 + b * math.cos(alphai_1)
            y1_i = y1 + b * math.sin(alphai_1)

            x2_i = x2 + b * math.cos(alphai_2)
            y2_i = y2 + b * math.sin(alphai_2)

            d = dist(x1_i, y1_i, x2_i, y2_i)
            Nleft = N - 2 * (i+1)
            dr_left = Nleft * b
            if d < dr_left:
                print "normal step, dist", d, " dist left = ", dr_left, " 1:x,y ", x1_i, y1_i, " 2:x,y ", x2_i, y2_i
                x1, y1 = x1_i, y1_i
                x2, y2 = x2_i, y2_i
            elif np.allclose(d, dr_left):
                # connect with a straight line
                print "straight line step, dist", d, " dist left = ", dr_left
                x1, y1 = x1_i, y1_i
                x2, y2 = x2_i, y2_i
                for j in range(Nleft):
                    x1_i = x1 + b * i
                    y1_i = y1 + b * i
                    X[i+j+1], Y[i+j+1] = x1_i, y1_i
                    print "steps ", x1_i, y1_i
                break
            elif d > dr_left:
                # do some magic with reiterating the last step
                print "too close line step, dist", d, " dist left = ", dr_left
                # x1, y1 = x1_i, y1_i
                # x2, y2 = x2_i, y2_i
                x = np.linspace(x1, x2, Nleft, endpoint=False, retstep=False)
                y = np.linspace(y1, y2, Nleft, endpoint=False, retstep=False)
                for j in range(Nleft):
                    X[i+j+1], Y[i+j+1] = x[j], y[j]
                    print "steps ", x[j], y[j]
                break
            else:
                print "error"


        print X
        print Y
        plt.plot(X, Y, 'bo-')


# for x1 in [0]:
#     for x2 in [2, 8]:
#         if x1<x2:
#             y1 = 0
#             y2 = x2

#             build_chain(x1, y1, x2, y2, b, N)

if __name__ == '__main__':
    x1 = np.array([2, 0, 0])
    x2 = np.array([6, 0, 2])
    bounds = 10*np.ones(3)
    b = 1.
    N = 40

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(x1[0],x1[1],x1[2], 'ro',color='red',marker='D',s=50)
    ax.scatter(x2[0],x2[1],x2[2], 'ro',color='red',marker='D',s=50)

    # (X, Y, Z) = build_chain_3d_final(x1, x2, N, b, bounds)
    Coord, _ = build_chain_3d_final(x1, x2, N, b, bounds)
    # ax.scatter(xx1, yy1, zz1, 'ro',color='red',marker='D',s=50)
    # ax.scatter(xx2, yy2, zz2, 'ro',color='red',marker='D',s=50)

    # ax.plot(X[:2], Y[:2], Z[:2], 'rD-')
    # ax.plot(X[-2:], Y[-2:], Z[-2:], 'rD-')
    # ax.plot(X[1:-1], Y[1:-1], Z[1:-1], 'bo-')
    # ax.plot(Coord[1:-1,0], Coord[1:-1,1], Coord[1:-1,2], 'bo-')
    ax.plot(Coord[:,0], Coord[:,1], Coord[:,2], 'bo-')

    plt.show()


