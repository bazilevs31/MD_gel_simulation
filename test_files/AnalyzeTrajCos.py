#!/usr/bin/env python

import numpy as np
from MDAnalysis import *
import os
from read_parameters import read_traj_vmd
from numba import double, long_, int_, jit, autojit
# from AnalyzeChain import get_bondlist_coords
from save_plots import save_plot_cryst_cos

#Trajectory analysis
# 1) analyze g2 alignment parameter 
# 2) analyze gyr parameter
# 3) analyze R2 end to end distance
# 4) analyze R2/gyr^2 parameter
# plot and store all variables to the figures folder
# 
def get_bonds(u):
    """given universe return chords of the (i+1-i-1) bonds and coords of central i atoms"""
    angles = u.angles
    chords = angles.atom3.positions - angles.atom1.positions 
    coords = angles.atom2.positions
    norm = np.linalg.norm(chords,axis=1)
    chords /= norm[:, None] #the norm vector is a (nx1) and we have to create dummy directions -> (n,3)
    return chords,coords
def get_angles(u):
    """
    given universe return chords of the (i+1-i-1) atoms(central ith atom) and angles
    input: universe
    output: atoms(coordinates), angles(angle_values)
    """
    angles = u.angles
    angle_values = angles.angles()
    atoms = angles.atom2.positions
    # atoms = np.ascontiguousarray(atoms, dtype=np.float32), 
    # angle_values = np.ascontiguousarray(angle_values, dtype=np.float32)
    atoms = np.asarray(atoms,dtype=np.float32)
    angle_values = np.asarray(angle_values,dtype=np.float32)
    return atoms, angle_values


# @jit('float32(float32[:,::1],float32[:],float32)',nopython=True)
@jit
def Calcg2(atoms, angles,L):
    """
    calculate g2 parameter
    input: atoms Natoms*3 array, angles Natoms*1 array of angles
    here Natoms is not number of atoms in the system, is number of angle vertices
    """
    # d local distance parameter , tmp needed to calculate d
    # dPhi - local angle difference, dist_around = cutoff
    # cos2 - local g2 parameter , n - counter of neighbours
    # g2 - the result of our program 

    Natoms = atoms.shape[0]
    M = atoms.shape[1]
    cos2 = 0.0; g2 = 0.0; d = 0.0; tmp = 0.0; dPhi = 0.0; dist_around = 8.0 ; n = 0
    for i in range(Natoms):
        n = 0
        cos2 = 0.0
        for j in range(Natoms):
            if i!=j:
                d = 0.0
                tmp = 0.0
                for k in range(M):
                    tmp = atoms[i,k] - atoms[j,k]
                    d += tmp*tmp
                d = np.sqrt(d)
                if d<dist_around:
                    dPhi = 2.0*(angles[i]-angles[j])
                    cos2 += np.cos(dPhi)
                    n += 1 #number of neighbours
        g2 += cos2/float(n)
        # print cos2/float(n)
    g2 /= Natoms
    return g2

def main():
    args = read_traj_vmd()
    psffile = os.path.splitext(args.psffile)[0]
    u = Universe(psffile+'.psf', args.traj)
    time, g2 = [], []

    box = u.universe.dimensions[:3]
    L = box[0]
    print (L)

    for ts in u.trajectory[args.startframe:args.endframe:args.trajskip]:
        # bonds, atoms = get_bondlist_coords(u)
        atoms,angles = get_angles(u)

        # bonds, atoms = get_bonds(u)
        N=atoms.shape[0]
        g = 0.0
        # print (atoms,angles)
        g = Calcg2(atoms,angles,L)
        # result =  ft.sparam(natoms=N,bonds=bonds,atoms=atoms,around=2.0,s=s)
        print "frame is " , ts.frame, " order = ", g
        time.append(ts.frame)
        g2.append(g)
        
        # Rgyr = np.sum(np.array([myres.atoms.radiusOfGyration(pbc=True) for myres in u.residues]))
        # Rgyr /= len(u.residues)
        # gyr.append(Rgyr*Rgyr)

        # R2 = np.sum(np.array([distance_sq_coords(myres.atoms[0].pos,myres.atoms[-1].pos) for myres in u.residues]))
        # R2 /= len(u.residues)
        # gr2.append(R2)

    time = np.array(time)
    g2 = np.array(g2)
    
    save_plot_cryst_cos(time,g2,psffile)

if __name__ == '__main__':
    main()


# print(os.getcwd()+"/" + "\n")
# to get the current directory
#
#npz=np.load('../coords.txt.npz')
#atoms=npz['atoms']
#bonds=npz['bonds']
