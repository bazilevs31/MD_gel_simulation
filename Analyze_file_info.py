#!/usr/bin/python

"""Summary
"""
import re
import pandas as pd
import numpy as np
import string
from gridData import Grid
import matplotlib.pyplot as plt
import os


def plot_axes_get_max_min(x, y, ax, **kwargs):
    """plots x and y with params"""
    ax.cla()
    ax.set(xlabel='{0[xlabel]}'.format(kwargs), ylabel='{0[ylabel]}'.format(kwargs),title='{0[title]}'.format(kwargs))
    # ax.titlesize==12

    # with plt.style.context(('seaborn')):
        # plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')
    # plt.show()
    ax.plot(x, y, kwargs['style'])
    print(kwargs['ylabel'])
    if 'phiz' in kwargs['ylabel']:
        print('phiz')
        mymax, mymin = y[3:-4].max(), y[3:-6].min()
    else:
        mymin, mymax = y[2:-2].min(), y[2:-3].max()

    # with plt.style.context(('seaborn')):
    ax.axhline(y=mymax,c='k',ls='-',lw=1.9)
    ax.axhline(y=mymin,c='k',ls='-',lw=1.9)

    return mymax, mymin


def get_min_max(_x, _y):
    """gets min max of the data"""
    Amp_1 = np.min(_y[1])
    Amp_2 = np.min(_y[1])
    return Amp_1, Amp_2


def get_steps_from_profile(_x, _y, ax, fig, variable='Plat', name='profile'):
    """tweaks data to make it fittable
    could be
    phi
    Plat
    concentration

    Args:
        _x (numpy array): the z coordinate of the array
        _y (numpy array): the parameter of interest
        ax (matplotlib axis): that I will use for plotting
        fig (matplotlib figure): that will be used for plotting
        variable (str, optional): the variable analyzed (depending on it the delta steps will be different)
        name (str, optional): the name of the profile

    Returns:
        float: mymax
        float: mymin
    """
    mymax, mymin = plot_axes_get_max_min(_x, _y, ax, style='.--', xlabel='z', ylabel=variable, title=name)
    directory='./pdf'
    if not os.path.exists(directory):
        os.makedirs(directory)
    fig.savefig('{directory}/{mytype}_{name}.png'.format(directory=directory, name=name, mytype=variable))
    # fig.clf()
    return mymax, mymin


def get_phi_method_2(_y, temp, n_init=4):
    """
    analyzes the density of the file using method 2
    """
    y_curve = temp * np.log(_y[n_init:]/_y[n_init])
    phi_method_2 = np.abs(np.min(y_curve) - np.max(y_curve))
    return phi_method_2


def get_ionization_from_string(_x):
    """
    analyzes the string and returns the ionization
    """
    if len(_x)>1:
        string_ionization = _x[0] + '.' + _x[1:]
        ionization = float(string_ionization)
    else:
        ionization = float(_x)
    return ionization

def get_temp_from_string(_x):
    """analyzes temp string and returns the float"""
    temperature = float(_x[0] + '.' +_x[1])
    return temperature


def get_profile_params(_filename):
    """reads profile file and analyzes which data file it belonged to"""
    # searchObj = re.search( r'(.*?)_sc(.*?)', _filename)
    print("analyzing profile name",_filename)
    searchObj = re.search(r'sc1844_nm100_l50_fl(.*?)_fr(.*?)_t(.*)', _filename)
    print(searchObj.groups())
    if searchObj:
        _fl = searchObj.group(1)
        _fr = searchObj.group(2)
        _temp = searchObj.group(3)
    fl = get_ionization_from_string(_fl)
    fr = get_ionization_from_string(_fr)
    print(_temp)
    temp = get_ionization_from_string(_temp)
    gel_name = 'sc1844_nm100_l50_fl{fl}_fr{fr}'.format(fl=_fl, fr=_fr)
    return fl, fr, temp, gel_name


def get_datname(_filename):
    """reads profile file and analyzes which data file it belonged to

    Args:
        _filename (string): filename with extension

    Returns:
        TYPE: Description
    """
    searchObj = re.search( r'(.*)_sc(.*?)_press.txt', _filename)
    if searchObj:
        found = searchObj.group(2)
    print(found)
    dataname = './data/sc'+found+'.dat'
    return dataname


def get_profile_filename(_filename):
    """reads profile file and analyzes which data file it belonged to

    Args:
        _filename (TYPE): Description

    Returns:
        TYPE: Description
    """
    # searchObj = re.search( r'(.*?)_sc(.*?)_ljeq1_coul_eq_coulsim_(.*?).txt', _filename)
    searchObj = re.search( r'(.*?)_sc(.*?)_ljeq1_coul_eq_coulsim_(.*?).txt', _filename)
    if searchObj:
        found1 = searchObj.group(1)
        found2 = searchObj.group(2)
        found3 = searchObj.group(3)
    print( found1, found2, found3)
    profiletype = found1
    # searchfilename=re.search(r'(.*)_ljeq1(.*)',found2)
    # dataname = 'sc'+searchfilename.group(1)
    # simtype = found3
    # searchfilename=re.search(r'(.*)_ljeq1(.*)',found2)
    dataname = 'sc'+found2
    simtype = found3
    print( "everything parsed  = ", profiletype, dataname, simtype)
    return profiletype, dataname, simtype

def get_info_dat(inFileName):
    """analyzes a .dat file in the ./data/....dat
    directory to get number of atoms and dimensions

    Args:
        inFileName (TYPE): Description

    Returns:
        TYPE: Description
    """
    inFile = open(inFileName, "r")
    lines = inFile.readlines()
    inFile.close()

    verbose = False

    for line in lines:
        if "atoms" in line:
            list_line = line.split(" ")
            N_atoms = int(list_line[0])
            if verbose:
                print( "N_atoms", N_atoms)
        elif "xlo xhi" in line:
            list_line = line.split(" ")
            xlo = float(list_line[0])
            xhi = float(list_line[1])
            if verbose:
                print( xlo, xhi)
        elif "ylo yhi" in line:
            list_line = line.split(" ")
            ylo = float(list_line[0])
            yhi = float(list_line[1])
            print( ylo, yhi)
        elif "zlo zhi" in line:
            list_line = line.split(" ")
            zlo = float(list_line[0])
            zhi = float(list_line[1])
            if verbose:
                print( zlo, zhi)
    return N_atoms, abs(xlo - xhi), abs(ylo - yhi), abs(zlo - zhi)




def analyze_pressure_profile(inFileName):
    """analyzes pressure profile and saves the output data onto
    csv file
    with press<name of atomtype><name pressure type>_<gel name>.csv

    Args:
        inFileName (string): name of the infput file with the extensions

    Returns:
        pandas dataframe: with the resulting information from the file
    """

    # inFileName = "Pressurealltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt"
    inFile = open(inFileName, "r")
    lines = inFile.readlines()
    inFile.close()


    # sc1844_nm100_l50_fl01_fr05_t11_ljeq1_coul_eq_coulsim.dat

    # find slab thickness (delta):
    for line in lines:
        if line[0] != '#': # ignore comments
            # words = string.split(line)
            words = line.split()
            if len(words) == 2 or len(words) == 3:
                nBins = int(words[1])

    # datname = get_datname(inFileName)
    # N_atoms, lx, ly, lz = get_info_dat(datname)
    # slabVolume = lx * ly * lz / (nBins - 1)
    # delta = lz / float(nBins - 1)

    x = np.zeros(nBins, dtype=[('Id', np.float32), ('Coord1', np.float32), ('Ncount', np.float32),  ('NumDensity', np.float32), ('Pressure', np.float32),('Pzz', np.float32),('Pxxyy', np.float32),('Pxx', np.float32),('Pyy', np.float32)])


    flag_reading = True
    traj_pd = []
    count = 0


    Lx = 200
    Ly = 200
    Lz = 950+50
    # volume = 218.*44.*44.
    volume = Lx*Ly*Lz
    bin_volume = volume / float(nBins)

    # calc & output P tensor components (file) and P_L-P_N (screen):
    # outFile = open('zP_xx_yy_zz', 'w')
    for line in lines:
        flag_reading = True
        if line[0] != '#': # ignore comments
            # words = string.split(line)
            words = line.split()
            if flag_reading and (len(words) != 3):
                # print( "reading Pressure data, count", count)
                x[count]['Id'] = int(words[0])
                x[count]['Coord1'] = float(words[1])
                x[count]['NumDensity'] = float(words[3])
                nCount = float(words[2])
                x[count]['Ncount'] = nCount
                cp1 = float(words[4])
                cp2 = float(words[5])
                cp3 = float(words[6])
                x[count]['Pressure'] = - nCount * (cp1 + cp2 + cp3) / (3 * bin_volume)
                x[count]['Pzz'] = - nCount * (cp3) / (bin_volume)
                x[count]['Pxx'] = - nCount * (cp1) / (bin_volume)
                x[count]['Pyy'] = - nCount * (cp2) / (bin_volume)
                x[count]['Pxxyy'] = - nCount * (cp1 + cp2) / ( 2*bin_volume)
                df = pd.DataFrame.from_records(x)
                # print( df)
                count += 1
            if len(words) == 3:
                # flag_reading = False
                # print( "len is 3")
                if count > 0:
                    traj_pd.append(df)
                    # traj_pandas = []
                    # flag_reading = False
                    frame = int(words[0])
                    count = 0

    profiletype, outputfile, simtype = get_profile_filename(inFileName)
    Combine_PD = pd.concat(traj_pd)
    Combine_PD.to_csv(profiletype+"_"+outputfile + '.csv', encoding='utf-8', index=False)
    # Combine_PD.to_csv(profiletype+"_"+outputfile + "_" + simtype+ '.csv', encoding='utf-8', index=False)
    return Combine_PD

def analyze_dx(_f):
    """analyzes dx _file

    Args:
        _f (string): filename with extension of the dx file that we want to analyze

    Returns:
        pandas: pandas dataframe columns=['coord_phi','phiz']
    """
    g = Grid(_f)
    nx, ny, nz = g.grid.shape
    Lx = g.edges[0][0]
    Ly = g.edges[1][0]
    Lzl = g.edges[2][0]
    Lzr = g.edges[2][-1]

    g_all = np.zeros(nz)

    count = 0

    for i in range(nx):
        for j in range(ny):
            g_all += g.grid[i][j]
            count += 1
    r_array = np.linspace(Lzl, Lzr, nz)

    yy = g_all/count

    A = np.vstack([r_array, yy])
    df = pd.DataFrame(A.T, columns=['coord_phi','phiz'])
    return df

# if __name__ == '__main__':
#     main()
