#!/usr/bin/env python

# Program: CreteGel.py
# Purpose: creates a polyelectrolyte gel
# Author:  Triandafilidi Vasiliy , PhD student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreateGel.py -h for help,

# Copyright (c) 2016 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version


import sys
import argparse
import numpy as np
import distances
import create_saw_vonmises_chain
class Lattice:
    """lattice class
    methods:
    gen_lattice
    cell
    property method:
    nsites

    input: ltype='sc', scale=1.0, n_=[1, 1, 1]
    ltype: fcc, diamond, bcc, sc
    scale: float unitcell length
    n_: array of [nx, ny, nz] dimensions

    to run : Lat = Lattice(), Lat.gen_lattice()
    output:
    self.box of dimensions
    number of sites can be accessed through the property Lat.nsites
    self.rcut : the first neighbor cutoff
    """
    def __init__(self, ltype='sc', scale=1.0, n_=[1, 1, 1]):
        self.ltype = ltype
        self.scale = scale
        self.n_ = n_
        self.c, self.rcut = self.cell(ltype, scale)
        self.site = []
        self.box = scale * (np.array(n_))

        # self.gen_lattice()

    def gen_lattice(self):
        site = disp = [0.0, 0.0, 0.0]
        for i in range(self.n_[0]):
            disp[0] = i * self.scale
            for j in range(self.n_[1]):
                disp[1] = j * self.scale
                for k in range(self.n_[2]):
                    disp[2] = k * self.scale
                    for a in self.c:
                        site = [a[l] + disp[l] for l in range(3)]
                        self.site.append(site)
        return self.box, self.rcut

    def __str__(self):
        sites = '\n'.join(map(lambda x: '%s: %s' % x, enumerate(self.site)))
        basis = '\n'.join(map(lambda x: '%s: %s' % x, enumerate(self.c)))
        box = '\n'.join(map(lambda x: '%s: %s' % x, enumerate(self.box)))
        doc_string = """\
        "Lattice with
        {nsites}
         type
         {ltype}
         basis is
         {basis}
         box =
         {box}
         sites =
         {sites}
         """.format(nsites=self.nsites, ltype=self.ltype, basis=basis, box=box, sites=sites)
        return doc_string
    @property
    def nsites(self):
        return len(self.site)

    @staticmethod
    def cell(ltype = 'fcc', scale = 1.0):
        """
        ltype: fcc, diamond, bcc, sc
        scale: float
        generate a unit cell of a given ltype
        returns basis
        """
        if ltype == 'diamond':
            basis = [[0.0, 0.0, 0.0],
                          [0.0, 0.5, 0.5],
                          [0.5, 0.0, 0.5],
                          [0.5, 0.5, 0.0],
                          [0.25, 0.25, 0.25],
                          [0.25, 0.75, 0.75],
                          [0.75, 0.25, 0.75],
                          [0.75, 0.75, 0.25]]
            # offset = 0.25
            offset = 0.0
            rcut = np.sqrt(3)/4.
        elif ltype == 'fcc':
            basis = [[0.0, 0.0, 0.0],
                      [0.0, 0.5, 0.5],
                      [0.5, 0.0, 0.5],
                      [0.5, 0.5, 0.0]]
            offset = 0.25
            rcut = np.sqrt(2)/2.
        elif ltype == 'bcc':
            basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.5]]
            offset = 0.25
            rcut = np.sqrt(2)/4.
        elif ltype == 'sc':
            basis = [[0.0, 0.0, 0.0]]
            offset = 0.5
            rcut = 1.
        else:
            print('unknown lattice, choices are: fcc, bcc, sc, diamond')
            sys.exit(1)
        # just shifting the origin
        for site in basis:
            for i in range(3):
                site[i] = (site[i] + offset) * scale
        return basis, rcut*scale



class Atom(object):
    """atom in a molecule or in a force field"""

    def __init__(self, name, m=1.0):
        self.name = name
        self.m = m
        self.neighbors = []
        self.q = 0
        self.ityp = 1
        self.chainnum = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.shift_idx = 0 # the resulting index will be idx + shift_idx
        self.idx = 0
        self._species = [] # to count what species have been created

    def __str__(self):
        # if hasattr(self, 'type'):
        #     return "atom {0:5s} {1:3s}  m = {2:7.3f}  q = {3:+7.4f}  "\
        #            "{4} {5}".format(self.name, self.type, self.m, self.q,
        #                           self.pot, str(self.par))
        # else:
        return "idx {idx} atom {0:5s}  m = {1:7.3f}  q = {2:+7.4f} neighs = {3:s} pos(x,y,z) {4:3.2f} {5:3.2f} {6:3.2f}\n".format(self.name, self.m, self.q, str(self.neighbors),self.x, self.y, self.z,idx=self.idx)

    @property
    def position(self):
        """coordinates of the atom
        Get the current cartesian coordinates of the atom.
        :Returns: a (3,) shape numpy array
        """
        # return np.array([self.x, self.y, self.z]).copy()
        return np.array([self.x, self.y, self.z])
        # return self.universe.coord.positions[self.index].copy()
    @position.setter
    def position(self, pos_xyz):
        """gets the numpy array of positions and sets the x y z
        """
        self.x = pos_xyz[0]
        self.y = pos_xyz[1]
        self.z = pos_xyz[2]
        return None
    @property
    def nneighs(self):
        return len(self.neighbors)


class Bond(object):
    """covalent bond in a molecule or in a force field"""

    def __init__(self, i = -1, j = -1, r = 0.0):
        self.i = i
        self.j = j
        self.r = r
        self.ityp = 1

    def __str__(self):
        if hasattr(self, 'name'):
            if self.i != -1:
                return "bond {0:5d} {1:5d}  {2}  {3} {4}".format(self.i + 1,
                        self.j + 1, self.name, self.pot, str(self.par))
            else:
                return "bond {0}  {1} {2}".format(self.name, self.pot,
                                                  str(self.par))
        else:
            # return "bond {0:5d}   ---{1:5d}   {2:1.4f}".format(self.i + 1, self.j + 1, self.r)
            return "bond {0:5d}   ---{1:5d}   {2}".format(self.i + 1, self.j + 1, self.r)



class Mol():
    """molecule"""

    def __init__(self,topol='diamond',  n_ = None, npoly=2, scale=1., ionization_every=1, csalt=0, densityfactor=2., valency=1):
        """
        plus additional parameter like chainnum which will distinguish different polymer chains as different molecules
        nodes and counter ions are 0
        the rest (polymer chains) are i = 1, 2 ...
        """
        self.nodes = []  # the nodes in the lattice <- class atom
        self.mons = []  # the monomers between the nodes <- class atom
        self.cions = [] # number of cions
        self.bonds = []

        self.neighs = [] # between nodes class bond
        self.npoly = npoly # length of the poly chain between nodes
        self.nnodes = 0
        self.cions = 0
        # self.nmons = 0
        self.topol = topol
        self.scale = scale
        self.box = [None]
        self.valency = valency # the valency of the backbone atoms
        # self._cell_dims
        self.name = 'gel'
        self.n_ = n_
        self.density_factor = float(densityfactor) # how different are the chain lengths
         # only every this atom will be charged
        # self.ionization_every = ionization_every
        self.ionization_every = ionization_every
        self.csalt = csalt
        self.nsalt = 0 #this will change
        self.verbose = False
        print "------  creating nodes ------- "
        self.create_nodes()
        print "------  creating mons ------- "
        self.create_mons()
        print "------  creating cions ------- "
        self.create_cions()
        if csalt > 0:
            print "------  creating salt macroions ------- "
            self.create_salt()

    def __str__(self):
        neighs = '\n'.join(map(lambda x: 'neigh: %s: %s-%s %s' % (x[0],x[1].i, x[1].j, x[1].r), enumerate(self.neighs)))
        bonds = '\n'.join(map(lambda x: 'bond: %s: %s-%s %s' % (x[0],x[1].i, x[1].j, x[1].r), enumerate(self.bonds)))
        header_string= '\nmolecule %s  %d atoms' % \
            (self.name, self.nnodes + self.nmons + self.ncions)
        parts_charges = '\n ionization alpha, frequency charges/parts = {0}/{1}\n'.format(self.ncharges, self.nmons + self.nnodes)
        parts_total_charge = '\n total charge of all particles {0}\n'.format(self.total_charge)
        parts_salt = '\nnumber salt molecules = {0}\n'.format(self.nsalt)
        parts_all = 'System, nodes = {nodes}, mons = {mons}, neighs = {neighs}, ncions = {ncions}, tot_charge = {tot_charge}'.format(nodes=self.nnodes, mons=self.nmons, neighs=self.nneighs, ncions=self.ncions, tot_charge=self.total_charge)
        return neighs + '\n' +"\n .... \n" +  bonds[:1000] + "\n .... \n" + header_string + parts_charges + parts_salt+parts_total_charge+parts_all

    @property
    def charge(self):
        q = 0.0
        for at in self.nodes:
            q += at.q
        for at in self.mons:
            q += at.q
        return q

    @property
    def total_charge(self):
        """
        use only in the end
        """
        q = 0.0
        for at in self.nodes:
            q += at.q
        for at in self.mons:
            q += at.q
        for at in self.cions:
            q += at.q
        return q

    @property
    def nmons(self):
        return len(self.mons)
    @property
    def nparticles(self):
        return self.nmons + self.nnodes + self.ncions + 2 * self.nsalt
    @property
    def ncharges(self):
        total_num_charges = 0
        for ns in self.nodes:
            if (abs(ns.q) > 0):
                total_num_charges += 1
        for ms in self.mons:
            if (abs(ms.q) > 0):
                total_num_charges += 1
        return total_num_charges
    @property
    def nbonds(self):
        return len(self.bonds)
    @property
    def nneighs(self):
        return len(self.neighs)
    @property
    def cell_dims(self):
        # return self._cell_dims
        return self.box
    def create_nodes(self):
        """creates nodes atoms"""
        Lat = Lattice(self.topol, self.scale, self.n_)
        self.box, self.rcut=Lat.gen_lattice()
        if self.verbose:
            print "generated lattice with bounds", self.box
        # self._cell_dims, self.rcut = Lat.gen_lattice()
        self.nnodes = Lat.nsites
        if self.topol is 'sc':
            self.max_neighbors = 5
        else:
            self.max_neighbors = 5
        if self.verbose:
            print self.nnodes
        self.nodes = [None] * self.nnodes
        for i, site in enumerate(Lat.site):
            # self.nodes[i] = Atom("C" + str(i+1), m=1.0)
            self.nodes[i] = Atom("C", m=1.0)
            self.nodes[i].q = -1*(abs(self.valency))
            self.nodes[i].position = np.array(site)
            self.nodes[i].idx = i + 1

        for i1 in xrange(0, self.nnodes):
            for i2 in xrange(i1 + 1, self.nnodes):
                d_pbc = distances.dist_w_pbc(self.nodes[i1].position, self.nodes[i2].position, bounds=self.box)
                # d_pbc = distances.dist_wo_pbc(self.nodes[i1].position, self.nodes[i2].position, bounds=self.box)
                free_spots_flag = (self.nodes[i1].nneighs <= self.max_neighbors)
                if np.allclose(d_pbc, self.rcut) and i1 != i2 and free_spots_flag:
                    self.nodes[i1].neighbors.append(i2)  #iden-fy neigh nodes
                    self.nodes[i2].neighbors.append(i1)
                    self.neighs.append(Bond(i1, i2, d_pbc))
            if self.verbose:
                print "----------done node------------"
                print "node ", self.nodes[i1]

            # print
        return None

    def create_mons(self):
        # """creates monomers along an axis"""
        _curr_id = self.nnodes  # the last position of atoms (nodes)
        # its also the first position of mons
        ionization_count = 0
        chaincount = 0
        for neigh in self.neighs:
            chaincount += 1
            a1 = self.nodes[neigh.i]
            a2 = self.nodes[neigh.j]
            # Lbox = self.box[0]
            if self.verbose:
                print self.box
                print self.box[0]
            if a1.position[0] < self.box[0] / 2.:
                npolychain = int(self.density_factor * self.npoly)
                mons_list = [None]*(npolychain - 1)
            elif a1.position[0] >= self.box[0] / 2.:
                npolychain = self.npoly
                mons_list = [None]*(npolychain - 1)
            # xyz_vec, bond_length = distances.mons_alng_line(a1.position, a2.position, npolychain , self.box)
            # xyz_vec, bond_length = create_chain_pe.build_chain_3d_final(a1.position, a2.position, npolychain + 1, 1., self.box, self.verbose)
            xyz_vec, bond_length = create_saw_vonmises_chain.create_chain(a1.position, a2.position, npolychain + 1,  self.box)
            if self.verbose:
            # if True:
                print "size chain ", xyz_vec.shape, xyz_vec
            ionization_count += 1
            if (ionization_count % self.ionization_every == 0):
                a1.q = -1*(abs(self.valency))

            for i, xyz in enumerate(xyz_vec):
                _curr_id += 1
                mons_list[i] = Atom("O", m=1.0)
                ionization_count += 1
                if (ionization_count % self.ionization_every == 0):
                    mons_list[i].q = -1*(abs(self.valency))
                mons_list[i].position = xyz
                mons_list[i].idx = _curr_id
                mons_list[i].ityp = 2
                mons_list[i].chainnum = chaincount
                if i < (npolychain - 2):
                    self.bonds.append(Bond(_curr_id, _curr_id + 1, bond_length))
            if (ionization_count % self.ionization_every == 0):
                mons_list[i].q = -1*(abs(self.valency))
            ionization_count += 1
            # print ionization_count

            self.bonds.append(Bond(a1.idx, mons_list[0].idx, bond_length))
            self.bonds.append(Bond(a2.idx, mons_list[-1].idx, bond_length))
            self.mons.extend(mons_list)


        return None
    def create_cions(self):
        """
        generate the cions to have the total charge 0
        make sure you run this after you generate the gel
        # we will choose 1/6th of the atoms from the right side
        pick random ones and move them to the left
        """
        shift_idx = self.nnodes + self.nmons
        # self.ncions = self.nnodes + self.nmons
        # self.ncions = self.ncharges
        self.ncions = int(abs(self.charge)) # because now we have a valency
        self.cions = [None] * self.ncions
        _N = self.ncions
        _L = self.box[0]
        R = np.random.uniform(0, 1, (_N, 3))  * self.box

        if self.density_factor > 1.:
            right_indices = np.where(R[:,0] > _L/2.)[0]
            if self.verbose:
                print "box", self.box
                print "x side", _L
                print "total number of cions", _N, self.ncions
                print "fraction of the right indices ", len(right_indices)
                print "applying procedure"
                print R
                print right_indices
            # (1+x)*n_per_part = N_tot
            # n_per_part = N_tot / (1+x)
            # so n per one :
            # xN_tot/(1+x) + N_tot/(1+x)
            # N_tot/(1+x) - y = xN_tot/(1+x)
            # (N/2) * y =  N_tot/(1+x)
            # y = 2/(1+x)
            np.random.shuffle(right_indices)
            _fraction = (_N/2.)*(self.density_factor - 1.)/(self.density_factor+1)
            # R[right_indices[:int(_N/6)],0] -= _L/2.
            R[right_indices[:int(_fraction)],0] -= _L/2.
            right_indices = np.where(R[:,0] > _L/2.)[0]
            if self.verbose:
                print "fraction of the right indices ", len(right_indices)
        for i, x_i in enumerate(R):
            self.cions[i] = Atom("H")
            self.cions[i].q = 1.
            self.cions[i].position = x_i
            self.cions[i].idx = shift_idx + i + 1
            self.cions[i].ityp = 3



    @staticmethod
    def n_from_c_salt(csalt, l):
        """gets number particles required in this volume for this concentration"""
        v = l**3
        factor = 2.5725*0.01
        nsalt = int(v*csalt/factor)
        print "nsalt = ", nsalt
        return nsalt

    def create_salt(self):
        """
        generate the salt macro ions
        """
        self.nsalt = Mol.n_from_c_salt(self.csalt, self.box[0])
        if self.nsalt > 0:
            shift_idx = self.ncions + self.nnodes + self.nmons

            # self.ncions = self.nnodes + self.nmons
            # self.csalt = 10 # total is 200 : cations(+) and anions(-))
            self.cations = [None] * (self.nsalt)
            self.anions = [None] * (self.nsalt)
            R = np.random.uniform(0, 1, (self.nsalt, 3)) * self.box
            for i, x_i in enumerate(R):
                self.cations[i] = Atom("Na")
                self.cations[i].q = 1.
                self.cations[i].position = x_i
                self.cations[i].idx = shift_idx + i + 1
                self.cations[i].ityp = 4
            shift_idx += self.nsalt
            # R = np.random.uniform(0, 1, (self.nsalt, 3)) * self.box
            R = np.random.rand(self.nsalt, 3) * self.box
            for i, x_i in enumerate(R):
                self.anions[i] = Atom("Cl")
                self.anions[i].q = -1.
                self.anions[i].position = x_i
                self.anions[i].idx = shift_idx + i + 1
                self.anions[i].ityp = 5




def writexyz(mol, filename='mol_default.xyz'):
        outfile = (filename).rsplit('.', 1)[0] + '.xyz'
        nparticles = mol.nparticles
        # nbonds = mol.nbonds
        # for sp_pos in species_positions:
        #     nparticles += sp_pos.nparticles
        #     nbonds += sp_pos.nbonds
        with open(outfile, 'w') as f:
            # f.write(str(len(self.nodes)) + '\n')
            f.write(str(nparticles) + '\n')
            f.write('elementary file' + '\n')
            for a in mol.nodes:
                f.write("{0:5s} {1:15.6f} {2:15.6f} {3:15.6f}\n".format(\
                        a.name, a.x, a.y, a.z))
            for a in mol.mons:
                f.write("{0:5s} {1:15.6f} {2:15.6f} {3:15.6f}\n".format(\
                        a.name, a.x, a.y, a.z))
            for a in mol.cions:
                f.write("{0:5s} {1:15.6f} {2:15.6f} {3:15.6f}\n".format(\
                        a.name, a.x, a.y, a.z))
            if mol.csalt:
                for a in mol.cations:
                    f.write("{0:5s} {1:15.6f} {2:15.6f} {3:15.6f}\n".format(\
                            a.name, a.x, a.y, a.z))
                for a in mol.anions:
                    f.write("{0:5s} {1:15.6f} {2:15.6f} {3:15.6f}\n".format(\
                            a.name, a.x, a.y, a.z))
        return None

def write_data(mol, filename):
    """writes lammps data file"""
    outfile = (filename).rsplit('.', 1)[0] + '.data'
    nparticles = mol.nparticles
    nbond = mol.nbonds
    # xlow, xhigh = mol.cell_dims
    if mol.csalt > 0:
        ntypes = 5
    else:
        ntypes = 3

    header_string = """\
LAMMPS data file. CGCMM style. atom_style full generated by VMD/TopoTools v1.5 on Mon Aug 01 15:24:50 PDT 2016
{natoms} atoms
{nbonds} bonds
0 angles
0 dihedrals
0 impropers
{ntypes} atom types
1 bond types
0 angle types
0 dihedral types
0 improper types
{xlo} {xhi}  xlo xhi
{ylo} {yhi}  ylo yhi
{zlo} {zhi}  zlo zhi

Masses

1 {mass_nodes} # nodes
2 {mass_mons} # monomers
3 {mass_cions} # cions
""".format(natoms=nparticles, nbonds=nbond,ntypes=ntypes, xlo=0., ylo=0., zlo=0., xhi=mol.cell_dims[0], yhi=mol.cell_dims[1], zhi=mol.cell_dims[2], mass_nodes=mol.nodes[0].m, mass_mons=mol.nodes[0].m, mass_cions=mol.cions[0].m)
    if mol.csalt > 0:
        header_string += """\
4 {mass_cations} # cationNa+
5 {mass_anions} # anionCl-
    """.format(mass_cations=mol.cations[0].m, mass_anions=mol.anions[0].m)
    header_string += """\nAtoms\n"""


    nmol=0
    with file(outfile, 'w') as f:
        f.write(header_string + "\n")
        for ia, a in enumerate(mol.nodes):
            # print a
            f.write("{index} {nmol} {ityp} {charge} {x} {y} {z}\n".format(
                    index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                    charge=a.q, x=a.x, y=a.y, z=a.z))

        for ia, a in enumerate(mol.mons):
            # print a
            f.write("{index} {nmol} {ityp} {charge} {x} {y} {z}\n".format(
                    index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                    charge=a.q, x=a.x, y=a.y, z=a.z))
        for ia, a in enumerate(mol.cions):
            # print a
            f.write("{index} {nmol} {ityp} {charge} {x} {y} {z}\n".format(
                    index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                    charge=a.q, x=a.x, y=a.y, z=a.z))
        if mol.csalt > 0:
            for ia, a in enumerate(mol.cations):
                # print a
                f.write("{index} {nmol} {ityp} {charge} {x} {y} {z}\n".format(
                        index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                        charge=a.q, x=a.x, y=a.y, z=a.z))
            for ia, a in enumerate(mol.anions):
                # print a
                f.write("{index} {nmol} {ityp} {charge} {x} {y} {z}\n".format(
                        index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                        charge=a.q, x=a.x, y=a.y, z=a.z))


        f.write('\nBonds\n\n')

        for ib, b in enumerate(mol.bonds):
            ib += 1
            # print b
            # f.write(" {0:7d} {1:7d}".format(b.j + 1, b.i + 1))
            f.write("{bond_id} {bond_type} {j:7d} {i:7d} \n".format(bond_id=ib,bond_type=b.ityp, j=b.j, i=b.i))

def read_cmd_line():
    """parses the command line parameters and return the args"""
    parser = argparse.ArgumentParser(description='program to create nodes, counter ions and salt ions (both cations and anions)')
    parser.add_argument('-f', '--filename', default='mol')
    parser.add_argument('-lt', '--lattice', default='diamond', help='generate the lattice of this type')
    parser.add_argument('-nm', '--nmonomers', default=40, type=int)
    parser.add_argument('-val', '--valency', default=1, type=int,help='valency of the backbone atoms')
    parser.add_argument('-l', '--length_scale', default=100., type=float, help='for diamond lattice the neighbor is ~0.4*scale, so for bondlength 1.0 one needs 4 monomers for length_scale=10.')
    parser.add_argument('--auto_length', dest='auto_length', action='store_true', help='automatically calculate the chain length based on the equilibrium length of the Flory chain in a good solvent Rnu = bN^nu, nu = 0.588, for diamond lattice the neighbor is Rnu = l = N^0.588/0.4')
    parser.add_argument('-a', '--alpha', default=1., type=int, help='degree of ionization, only 1,alpha+1,2alpha+1.... will be charged, 1 - half of the monomers will be charged')
    parser.add_argument('-df', '--density_factor', default=1., type=float, help='degree of difference, by giving the chain length scaling factor, chains on the left are nm*densityfactor long on the right nm long')
    # parser.add_argument('-s', '--salt', default=0., type=float, help='concentration of salt ions, no macro-ions if 0.')
    parser.add_argument('-s', '--salt', default=0, type=float, help='salt concentration in [M=mole/litre] so 1 M = 0.025725 particles/sigma^3')
    # parser.add_argument('-d', '--dimensions', default=[3, 3, 3])
    parser.add_argument('-d','--dimensions', nargs='+', type=int, default=[1, 1,1])
    args = parser.parse_args()
    if args.auto_length == True:
        print "doing autolength"
        scale = float(args.density_factor)
        n = int(args.nmonomers)
        nscaled = n * scale
        l_flory = 1. * nscaled**0.588
        args.length_scale=l_flory/.4
        print "flory length = ", l_flory
        print "length set to flory_length/.4 (for diamond neigh distnace) = ", l_flory/.4
    else:
         print "nothing done automatically"
    print"Check that args.do_something=" + str(args.auto_length) + " is always a bool"
    return args

#vmd
def write_vmd(filename, length_scale):
    """writes the vmd file for vizualization"""
    input_data = (filename).rsplit('.', 1)[0] + '.data'
    viz_file = (filename).rsplit('.', 1)[0] + '.vmd'
    r_vdw = length_scale / 25.
    r_cpk = length_scale / 36.
    r_salt = length_scale / 72.
    r_point = 3.5
    header_string="""\
color Display Background white
topo readlammpsdata {input_data} full waitfor all
#nodes and mons
mol modselect 0 0 name is "1" or name is "2"
mol modcolor 0 0 ColorID 0
mol modstyle 0 0 VDW {r_cpk} 12.000000


#cions
mol addrep 0
mol modselect 1 0 name is "3"
mol modcolor 1 0 ColorID 7
mol modstyle 1 0 Points {r_point}

# # Na+ cation
# mol addrep 0
# mol modselect 2 0 name is "4"
# mol modcolor 2 0 ColorID 27
# mol modstyle 2 0 VDW {r_salt} 12.000000

# # Cl- anion
# mol addrep 0
# mol modselect 3 0 name is "5"
# mol modcolor 3 0 ColorID 16
# mol modstyle 3 0 VDW {r_salt} 12.000000
# VDW nodes

mol addrep 0
mol modselect 2 0 name is "1"
mol modcolor 2 0 ColorID 1
mol modstyle 2 0 VDW {r_vdw} 12.000000

set sel_all [atomselect top "name is 3"]
set minmax [measure minmax $sel_all]
set box [vecscale 1.0 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
pbc set $box
pbc box
display resetview
render TachyonInternal {input_data}.png /usr/bin/open %s

""".format(input_data=input_data, r_vdw=r_vdw, r_cpk=r_cpk, r_point=r_point, r_salt=r_salt)
    with file(viz_file, 'w') as f:
        f.write(header_string)

def main():
    """program to generate the  pe lattice that can be used as a module"""
    args = read_cmd_line()
    verbose = False
    print args
    Gel = Mol(topol=args.lattice, scale=args.length_scale, npoly=args.nmonomers, n_=args.dimensions, ionization_every=args.alpha, csalt=args.salt,densityfactor=args.density_factor, valency=args.valency)
    if verbose:
        print Gel
    print Gel
    writexyz(Gel, args.filename)
    write_data(Gel, args.filename)
    write_vmd(args.filename, args.length_scale)

    # print "System, nodes = {0}, mons = {1}, neighs = {2}".format(Gel.nnodes, Gel.nmons, Gel.nneighs)
if __name__ == '__main__':
    main()
