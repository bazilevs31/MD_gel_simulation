#!/usr/bin/env python

# Program: CreateDefChain.py
# Purpose: Creates monodisperse chains
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreateDefChain.py -h for help,
# example: python CreateDefChain.py -nch 100 -l 40 , for the system of 100 of C40
# Requires: read_parameters.py

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

from __future__ import print_function
# import os
import numpy as np
import read_parameters


def main():

    """
    write created distribution to a file
    input Nchains, ChainLength
    output ( def.chain)
    """
    args = read_parameters.read_create_mono()
    filename = 'def'
    parameters = {'randseed': np.random.randint(int(102220*args.Nchains)),
                'Nchains': args.Nchains,
                'ChainLength': args.ChainLength,
                'filename': filename
    }

    script = """\
Polymer chain definition

0.8442          rhostar
{0[randseed]}          random # seed (8 digits or less)
1               # of sets of chains (blank line + 6 values for each set)
0               molecule tag rule: 0 = by mol, 1 = from 1 end, 2 = from 2 ends

{0[Nchains]}             number of chains
{0[ChainLength]}             monomers/chain
1               type of monomers (for output into LAMMPS file)
1               type of bonds (for output into LAMMPS file)
0.85            distance between monomers (in reduced units)
1.05            no distance less than this from site i-1 to i+1 (reduced unit)

    """.format(parameters)

    with open('{0[filename]}{0[ChainLength]}n{0[Nchains]}.chain'.format(parameters), 'w') as f:
        f.write(script+'\n')
    return None

if __name__ == '__main__':
    main()
