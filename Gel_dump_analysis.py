#!/usr/bin/env python

import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
import re
import os
import argparse

import pandas as pd

def EOF(f):
    current_pos = f.tell()
    file_size = os.fstat(f.fileno()).st_size
    return current_pos < file_size




def analyze_pressure_dump(filename, Lx=200., Ly=200, Lz=900., N=10, bin_divide_flag=False, Natoms=113579):
    """
    analyzes pressure of a given
    """
    myfile = open(filename+'.txt')
    trajectory = []
    traj_pd = []
    frames = []

    for _ in range(3):
        next(myfile)
    count = 0
    while EOF(myfile):
        count += 1
        s = next(myfile) # info with the time step

        x = np.zeros(N, dtype=[('Chunk',np.float32), ('Coord1',np.float32), ('Ncount',np.float32),  ('density',np.float32), ('temp',np.float32), ('vx',np.float32), ('fx',np.float32),('c_pciKE[1]',np.float32), ('c_pciKE[2]',np.float32), ('c_pciKE[3]',np.float32), ('c_pciVIR[1]',np.float32), ('c_pciVIR[2]',np.float32), ('c_pciVIR[3]',np.float32), ('c_pgelELAS[1]',np.float32), ('c_pgelELAS[2]',np.float32), ('c_pgelELAS[3]',np.float32), ('c_pgelVIR[1]', np.float32), ('c_pgelVIR[2]', np.float32), ('c_pgelVIR[3]', np.float32), ('c_pgelPAIR[1]', np.float32), ('c_pgelPAIR[2]', np.float32), ('c_pgelPAIR[3]', np.float32)])

# Chunk Coord1 Ncount density/number temp vx fx c_pciKE[1] c_pciKE[2] c_pciKE[3] c_pciVIR[1] c_pciVIR[2] c_pciVIR[3] c_pgelELAS[1] c_pgelELAS[2] c_pgelELAS[3] c_pgelVIR[1] c_pgelVIR[2] c_pgelVIR[3] c_pgelPAIR[1] c_pgelPAIR[2] c_pgelPAIR[3]

        list_line = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)
        frame, _, _ = list_line
        frames.append(int(frame))
        # print( "reading lines")

        for i in xrange(N):
            count += 1
            s = next(myfile)
            list_line = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)
            # print( "reading line", i, list_line)
            for il, l in enumerate(list_line):
                x[i][il] = float(l)

        trajectory.append(x)

        # names = x.dtype.fields.keys()
        # data = x.dtype.fields.values()

        df = pd.DataFrame.from_records(x)
        traj_pd.append(df)

    myfile.close()



    # # volume = 218.*44.*44.
    volume = Lx*Ly*Lz
    # N_atoms = 113579
    # if bin_divide_flag:
    #     bin_volume = volume / float(N)
    # else:
    #     bin_volume = 1.

    bin_volume = volume / float(N)
    # bin_volume = volume
    # bin_volume /= float(Natoms)

    Combine_PD = pd.concat(traj_pd)
    FINAL_PD = pd.DataFrame()

    FINAL_PD['Coord1'] = Combine_PD['Coord1']
    FINAL_PD['p_ciKE'] = -1 * Combine_PD['Ncount'] * (Combine_PD['c_pciKE[1]'] + Combine_PD['c_pciKE[2]'] + Combine_PD['c_pciKE[3]'])/(3.*bin_volume)
    FINAL_PD['p_ciVIR'] = -1 * Combine_PD['Ncount'] * (Combine_PD['c_pciVIR[1]'] + Combine_PD['c_pciVIR[2]'] + Combine_PD['c_pciVIR[3]'])/(3.*bin_volume)
    FINAL_PD['p_gelELAS'] = -1 * Combine_PD['Ncount'] * (Combine_PD['c_pgelELAS[1]'] + Combine_PD['c_pgelELAS[2]'] + Combine_PD['c_pgelELAS[3]'])/(3.*bin_volume)

    FINAL_PD['p_gelVIR'] = -1 * Combine_PD['Ncount'] * (Combine_PD['c_pgelVIR[1]'] + Combine_PD['c_pgelVIR[2]'] + Combine_PD['c_pgelVIR[3]'])/(3.*bin_volume)
    FINAL_PD['p_gelPAIR'] = -1 * Combine_PD['Ncount'] * (Combine_PD['c_pgelPAIR[1]'] + Combine_PD['c_pgelPAIR[2]'] + Combine_PD['c_pgelPAIR[3]'])/(3.*bin_volume)

    # So now I have to
    # P_bin = (sigma_per_atom_xx + ... + sigma_per_atom_zz)/(bin_volume*3)
    # *N_atoms_per_bin
    # N_atoms_per_bin = number_density*N_atoms


    df_concat = FINAL_PD

    by_row_index = df_concat.groupby(df_concat.index)
    df_means = by_row_index.mean()
    by_row_index_2 = df_concat.groupby(df_concat.index)
    df_stds = by_row_index_2.std()

    # print( df_means.head())
    # print( df_stds.head())
    return df_means, df_stds


def analyze_energy_dump(filename, N=10):
    """
    analyzes energy dump
    """
    x = np.zeros(N, dtype=[('Chunk',np.float32), ('Coord1',np.float32), ('Ncount',np.float32),  ('density',np.float32),('temp',np.float32),('c_enciVIR',np.float32),('c_engelELAS',np.float32), ('c_engelPAIR',np.float32)])

# Chunk Coord1 Ncount density/number temp c_enciVIR c_engelELAS c_engelPAIR


    trajectory = []
    traj_pd = []
    frames = []

    myfile = open(filename+'.txt')
    # print( myfile)
    count = 0
    for line in myfile:
        li=line.strip()
        if not li.startswith("#"):
            list_line = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)
            if len(list_line) > 3:
                print( "count = ", count)

                print( "reading line", list_line)
                for il, l in enumerate(list_line):
                    x[count][il] = float(l)
                count += 1
            elif len(list_line) == 3:
                count = 0
                frame, _, _ = list_line
                trajectory.append(x)
                df = pd.DataFrame.from_records(x)
                traj_pd.append(df)
                frames.append(int(frame))


    Combine_PD = pd.concat(traj_pd)
    FINAL_PD = Combine_PD

    #     new_traj_pd.append(new_pd)

    # FINAL_PD = pd.concat(new_traj_pd)


    # df_concat = pd.concat((df1, df2))

    df_concat = FINAL_PD

    by_row_index = df_concat.groupby(df_concat.index)
    df_means = by_row_index.mean()
    by_row_index_2 = df_concat.groupby(df_concat.index)
    df_stds = by_row_index_2.std()

    print( df_means.head())
    print( df_stds.head())
    return df_means, df_stds




def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def analyze_density_dump(filename, N=10):

    x = np.zeros(N, dtype=[('Chunk',np.float32),('OrigID',np.float32), ('Coord1',np.float32), ('Ncount',np.float32),  ('density',np.float32)])


# Chunk Coord1 Ncount density/number temp c_enciVIR c_engelELAS c_engelPAIR


    trajectory = []
    traj_pd = []
    frames = []

    myfile = open(filename+'.txt')
    # print( myfile)
    count = 0
    for line in myfile:
        li=line.strip()
        if not li.startswith("#"):
            list_line = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)
            if len(list_line) > 3:
                print( "count = ", count)

                print( "reading line", list_line)
                for il, l in enumerate(list_line):
                    x[count][il] = float(l)
                count += 1
            elif len(list_line) == 3:
                count = 0
                frame, _, _ = list_line
                trajectory.append(x)
                df = pd.DataFrame.from_records(x)
                traj_pd.append(df)
                frames.append(int(frame))


    Combine_PD = pd.concat(traj_pd)
    FINAL_PD = Combine_PD

    #     new_traj_pd.append(new_pd)

    # FINAL_PD = pd.concat(new_traj_pd)


    # df_concat = pd.concat((df1, df2))

    df_concat = FINAL_PD

    by_row_index = df_concat.groupby(df_concat.index)
    df_means = by_row_index.mean()
    by_row_index_2 = df_concat.groupby(df_concat.index)
    df_stds = by_row_index_2.std()

    print( df_means.head())
    print( df_stds.head())
    return df_means, df_stds

