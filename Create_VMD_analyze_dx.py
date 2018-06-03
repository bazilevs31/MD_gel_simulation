#!/usr/bin/env python

# Program: Create_VMD_analyze_dx.py
# Purpose: creates a vmd script to analyze the system
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreCreate_VMD_analyze_dx.py -h for help,

# Copyright (c) 2018 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import os
import fnmatch
import argparse




def write_file(datafile, dcdfile,  dxfile, molnumber, vmdfile, ewaldfactor=3., nx=8, ny=8, nz=18):
    """writes the vmd file"""

    headstring = """
    package require pbctools
    package require topotools"""

    params  = { 'datafile': datafile,
                'dcdfile' : dcdfile,
                'dxfile' : dxfile,
                'vmdfile' : vmdfile,
                'count': molnumber}

    body_string = """
    topo readlammpsdata {0[datafile]} full waitfor all

    animate read dcd {0[dcdfile]} waitfor all
    package require pmepot
    pbc wrap
    # # pmepot -mol 0  -ewaldfactor 0.25
    pmepot -mol {0[count]} -ewaldfactor {0[ewaldfactor]}  -frames all  -grid {{ {0[nx]} {0[ny]} {0[nz]} }}    -dxfile {0[dxfile]}.dx -loadmol none
    mol delete {0[count]}

    """.format(params)


    with open("{0[vmdfile]}.vmd".format(params), 'a') as f:
        f.write(headstring  + body_string +  '\n' +'\n')
    return None


def walk_and_create_vmd():
    """walks through the director ./data and analyzes the files"""
    dcdmatches = []
    datamatches = []
    dxmatches = []

    datadir_to_walk = './data/'
    press_filter = 'sc*press.dcd'
    dcd_suffix = '_ljeq1_coul_eq_coulsim_press'
    _vmdfile = 'glob_vmd_all'

    for root, dirnames, filenames in os.walk(datadir_to_walk):
        for filename in fnmatch.filter(filenames, press_filter):
            dcdmatches.append(os.path.join(root, filename))
            base=os.path.basename(filename)
            cutname = os.path.splitext(base)[0]
            _dxfile = cutname.replace(dcd_suffix,'')
            dxmatches.append(_dxfile)
            print( cutname)
    # print( dcdmatches)

    for _f in dxmatches:
        for root, dirnames, filenames in os.walk(datadir_to_walk):
            data_filter = '{filename}_{name_suffix}.data'.format(filename=_f, name_suffix='_ljeq1_coul_eq_coulsim')
            for filename in fnmatch.filter(filenames, data_filter):
                datamatches.append(os.path.join(root, filename))

    # print( datamatches)
    # print( dcdmatches)
    # print( dxmatches)

    with open("{vmdfile}.vmd".format(vmdfile=_vmdfile), 'w') as f:
        f.write('\n')

    for i,(_dat, _dcd, _dx) in enumerate(zip(datamatches, dcdmatches, dxmatches)):
        print( i)
        print( _dat)
        print( _dcd)
        print( _dx)
        write_file(_dat, _dcd, _dx, i, _vmdfile)

    with open("{vmdfile}.vmd".format(vmdfile=_vmdfile), 'a') as f:
        f.write(' \n exit \n')
    return None

if __name__ == '__main__':

    def is_valid_file(parser, arg):
        """
        Check if arg is a valid file
        """
        arg = os.path.abspath(arg)
        if not os.path.exists(arg):
            parser.error("The file %s doesn't exist " % arg)
        else:
            return arg


    # parser = ArgumentParser(description="""First paragraph

    #                                        Second paragraph

    #                                        Third paragraph""",
    #                                        usage='%(prog)s [OPTIONS]',
    #                                        formatter_class=RawTextHelpFormatter)

    # # parser.add_argument("--csv_with_files", dest="csv_with_files",
    # #                         default="new_gel_params.csv",
    # #                         type=lambda x: is_valid_file(parser, x),
    # #                         help="Input trajectory file)", metavar="FILE")

    # # parser.add_argument('-a', action='append', dest='files',
    # #                     type=str,
    # #                     default=[],
    # #                     help='Add values to analyze values to a list type (default: %(default)s)',
    #                     # )
    # parser.add_argument("-datas",
    #                     "--data_suffix",
    #                     dest="data_suffix",
    #                     default='_ljeq1_coul_eq_coulsim',
    #                     type=str,
    #                     help="suffix for the datafile gelname_<suffix> (default: %(default)s)")
    # options = parser.parse_args()

    walk_and_create_vmd()
