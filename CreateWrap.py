#!/usr/bin/env python
import numpy as np
from MDAnalysis import *
from read_parameters import read_traj_vmd,create_wrapdcd
import os
from termcolor import colored


def main():
    """creates wraped trajectory"""
    args = read_traj_vmd()
    # psffile = os.path.splitext(args.psffile)[0]
    create_wrapdcd(args)
if __name__ == '__main__':
    main()