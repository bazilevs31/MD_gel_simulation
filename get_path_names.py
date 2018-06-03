#!/usr/bin/env python

# Program: get_path_names.py
# Purpose: different cool features of python to work with paths
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  import get_path_names

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


import os


def get_filename(s):
    """
    input: filename, maybe with a path,
    returns only filename
    """
    (_, filename) = os.path.split(s)  # get rid of the path
    (filename, _) = os.path.splitext(filename)  # get rid of the ext.
    return filename


def get_path(s):
    """
    input: string that represents file, with a path
    returns: path
    """
    (path, _) = os.path.split(s)  # get rid of the filename
    return path


def combine_in_path(*parts):
    """
    gets a list with whole bunch of stuff to combine
    combines them and produces a file with full path
    """
    fullpath = os.path.join(*parts)
    fullpath = os.path.normpath(fullpath)
    return fullpath


def main():
    """
    main
    """
    print get_filename('/lalal/alald/filemy.psf')
    # # this splits the path and the file

    # os.path.split(path)

    # # this splits the file and the file extension
    # os.path.splitext(path)

    # os.path.join(*parts)

    # os.path.normpath(path)
    return None

if __name__ == '__main__':
    main()
