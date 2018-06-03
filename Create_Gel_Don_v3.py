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
    def __init__(self, ltype='sc', scale=1.0, n_=[1, 1, 1], init_coord=[0, 0, 0]):
        self.ltype = ltype
        self.scale = scale
        self.n_ = n_
        self.c, self.rcut = self.cell(ltype, scale)
        self.site = []
        self.box = scale * (np.array(n_))
        self.init_coord = init_coord

        # self.gen_lattice()

    def gen_lattice(self):
        site = disp = [0.0, 0.0, 0.0]
        # site = disp = self.init_coord
        for i in range(self.n_[0]):
            disp[0] = i * self.scale
            for j in range(self.n_[1]):
                disp[1] = j * self.scale
                for k in range(self.n_[2]):
                    disp[2] = k * self.scale
                    for a in self.c:
                        site = [a[l] + disp[l] + self.init_coord[l] for l in range(3)]
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
        self.boundary_type = '' # could be top, bottom, no
    def __str__(self):
        # if hasattr(self, 'type'):
        #     return "atom {0:5s} {1:3s}  m = {2:7.3f}  q = {3:+7.4f}  "\
        #            "{4} {5}".format(self.name, self.type, self.m, self.q,
        #                           self.pot, str(self.par))
        # else:
        return "idx {idx} atom {0:5s}  m = {1:7.3f}  q = {2:+7.4f} neighs = {3:s} pos(x,y,z) {4:3.2f} {5:3.2f} {6:3.2f} {interface}\n".format(self.name, self.m, self.q, str(self.neighbors),self.x, self.y, self.z,idx=self.idx,interface=self.boundary_type)

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

    def __init__(self,topol='diamond',  n_ = None, npoly=2, scale=1., ionization_every=1, csalt=0, densityfactor=2., valency=1, switch_xz=True, pbc_=[1,1,1], rand_walk_chain=False):
        """
        plus additional parameter like chainnum which will distinguish different polymer chains as different molecules
        nodes and counter ions are 0
        the rest (polymer chains) are i = 1, 2 ...
        """
        self.npoly = npoly # length of the poly chain between nodes
        self.topol = topol
        self.scale = scale
        self.valency = valency # the valency of the backbone atoms
        self.n_ = n_
        self.ionization_every = ionization_every #only every this fraction of atoms will be charged
        self.csalt = csalt

        if self.topol == 'diamond':
            self.max_connections = 4
        else:
            self.max_connections = 5
        self.slab_offset = 10
        self.nprobes = 0
        self.dummy_probe = True
        self.name = 'gel'
        self.verbose = False
        self.pbc_x = pbc_[0]==1  # if its 1 then True, else False
        self.pbc_y = pbc_[1]==1
        self.pbc_z = pbc_[2]==1
        self.rand_walk_chain = rand_walk_chain

        #this will change
        self.nodes = []  # the nodes in the lattice <- class atom
        self.mons = []  # the monomers between the nodes <- class atom
        self.cions = [] # number of cions
        self.bonds = []
        self.neighs = [] # between nodes class bond
        self.bottom_sites = []  # the atoms on the top and bottom
        self.top_sites = []
        self.box = [None]
        self.system_box = [None]
        self.name_species = []
        self.species = [] # this will have all the types of atoms we have
        self.nspecies = [] # number of atoms of the species
        self.mspecies = []
        # it will have self.top_sites, ...

        self.switch_xz = switch_xz
        self.nnodes = 0
        self.nmons = 0
        self.nsalt = 0
        self.ncions = 0
        self.nsalt = 0
        self.ionization_count = 0
        self.nprobes = 0
        print "------  creating nodes ------- "
        self.create_nodes()
        print "------  creating mons ------- "
        self.create_mons()
        print "------  creating cions ------- "
        # self.create_cions()
        if csalt > 0:
            print "------  creating salt macroions ------- "
            self.csalt = csalt
            self.create_salt()


        self.create_cions()

        if self.nprobes > 0:
            self.create_probe_atoms()
        elif self.dummy_probe:
            self.create_dummy_probe()
        else:
            print "no probe atoms added"

    def __str__(self):
        nodes_info = '\n'.join(map(lambda x: 'node %s: %s' % (x[0],x[1]), enumerate(self.nodes)))

        neighs = '\n'.join(map(lambda x: 'neigh: %s: %s-%s %s' % (x[0],x[1].i, x[1].j, x[1].r), enumerate(self.neighs)))
        bonds = '\n'.join(map(lambda x: 'bond: %s: %s-%s %s' % (x[0],x[1].i, x[1].j, x[1].r), enumerate(self.bonds)))
        header_string= '\nmolecule %s  %d atoms' % \
            (self.name, self.nparticles)
        parts_charges = '\n ionization alpha, frequency charges/parts = {0}/{1}\n'.format(self.ncharges, self.nparticles )
        parts_total_charge = '\n total charge of all particles {0}\n'.format(self.total_charge)
        parts_salt = '\nnumber salt molecules = {0}\n'.format(self.nsalt)
        parts_all = 'System, nodes = {nodes}, mons = {mons}, neighs = {neighs}, ncions = {ncions}, tot_charge = {tot_charge}'.format(nodes=self.nnodes, mons=self.nmons, neighs=self.nneighs, ncions=self.ncions, tot_charge=self.total_charge)
        return nodes_info+neighs + '\n' +"\n .... \n" +  bonds[:1000] + "\n .... \n" + header_string + parts_charges + parts_salt+parts_total_charge+parts_all

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
        Goes through the species in the molecule,
        goes through atoms of the species
        gets their charge and sums it up
        """
        q = 0.0
        for atom_type in self.species:
            for at in atom_type:
                q += at.q
        return q

    @property
    def nmons(self):
        return len(self.mons)
    @property
    def nparticles(self):
        """
        Goes through the species in the molecule,
        goes through atoms of the species
        gets their number and sums it up
        """
        nparticles = 0
        for n_atom_type in self.nspecies:
            if n_atom_type is not 0:
                nparticles += n_atom_type
        return nparticles
    @property
    def ncharges(self):
        """
        Goes through the species in the molecule,
        goes through atoms of the species
        gets their charge and if it is above zero it sums it up
        """
        total_num_charges = 0
        for atom_type in self.species:
            for atom in atom_type:
                if (abs(atom.q) > 0):
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
        shift_idx = self.nparticles
        Lat = Lattice(self.topol, self.scale, self.n_)
        self.box, self.rcut=Lat.gen_lattice()
        if self.verbose:
            print "generated lattice with bounds", self.box
        self.nnodes = Lat.nsites
        self.system_box = self.box
        if self.verbose:
            print self.nnodes
        self.nodes = [None] * self.nnodes
        _m_nodes = 1.
        # if we have nx ny nz lattice
        # first ny*nz sites are so from 0 to ny*nz - 1 - are bottom
        # the last ny*nz sites are top
        ninterface_nodes = int(self.n_[1]*self.n_[2])
        print "ninterface_nodes", ninterface_nodes

        for i, site in enumerate(Lat.site):
            # self.nodes[i] = Atom("C" + str(i+1), m=1.0)
            self.nodes[i] = Atom("C", m=_m_nodes)
            self.nodes[i].position = np.array(site)
            self.nodes[i].idx = shift_idx + i + 1
            if i <= (ninterface_nodes - 1):
                self.nodes[i].boundary_type = 'bottom'
            elif i >= (self.nnodes - ninterface_nodes):
                self.nodes[i].boundary_type = 'top'
            else:
                self.nodes[i].boundary_type = 'no'
            if (self.ionization_count % self.ionization_every == 0):
                self.nodes[i].q = -1*(abs(self.valency))
                self.ionization_count += 1


        for i1 in xrange(0, self.nnodes):
            for i2 in xrange(i1 + 1, self.nnodes):
                d_w_pbc = distances.dist_w_pbc(self.nodes[i1].position, self.nodes[i2].position, bounds=self.box)
                d_wo_pbc = distances.dist_wo_pbc(self.nodes[i1].position, self.nodes[i2].position, bounds=self.box)
                i2_pos_pbc = distances.get_pbc(self.nodes[i1].position, self.nodes[i2].position, bounds=self.box) # returns position of i2 with pbc
                free_spots_flag = (self.nodes[i1].nneighs <= self.max_connections)
                close_neigh_flag = np.allclose(d_wo_pbc, self.rcut)
                pbc_neigh_flag = np.allclose(d_w_pbc, self.rcut)
                connect_flag = False
                xlo = 0
                xhi = self.box[0]
                ylo = 0
                yhi = self.box[1]
                zlo = 0
                zhi = self.box[2]
                # print self.nodes[i2].position
                if close_neigh_flag and i1 != i2 and free_spots_flag:
                    connect_flag = True

                if pbc_neigh_flag and i1 != i2 and free_spots_flag and self.pbc_x:
                    if ((i2_pos_pbc[0] < xlo) or (i2_pos_pbc[0] > xhi)):
                        connect_flag = True

                if pbc_neigh_flag and i1 != i2 and free_spots_flag and self.pbc_y:
                    if ((i2_pos_pbc[1] < ylo) or (i2_pos_pbc[1] > yhi)):
                        connect_flag = True
                if pbc_neigh_flag and i1 != i2 and free_spots_flag and self.pbc_z:
                    if ((i2_pos_pbc[2] < zlo) or (i2_pos_pbc[2] > zhi)):
                        connect_flag = True

                if connect_flag:
                    self.nodes[i1].neighbors.append(i2)  #iden-fy neigh nodes
                    self.nodes[i2].neighbors.append(i1)
                    self.neighs.append(Bond(i1, i2, d_w_pbc))
            if self.verbose:
                print "----------done node------------"
                print "node ", self.nodes[i1]
        self.name_species.append('nodes')
        self.species.append(self.nodes)
        self.nspecies.append(self.nnodes)
        self.mspecies.append(float(_m_nodes))
        return None

    def create_mons(self):
        """creates monomers along an axis"""
        shift_idx = self.nparticles
        _curr_id = shift_idx  # the last position of atoms (nodes)
        chaincount = 0
        _m_mons = 1.
        for neigh in self.neighs:

            chaincount += 1
            a1 = self.nodes[neigh.i]
            a2 = self.nodes[neigh.j]
            npolychain = self.npoly
            mons_list = [None]*(npolychain - 1)

            if self.rand_walk_chain:
                xyz_vec, bond_length = distances.chain_alng_line(a1.position, a2.position, npolychain , self.box)
            else:
                xyz_vec, bond_length = distances.mons_alng_line(a1.position, a2.position, npolychain , self.box)
            self.ionization_count += 1

            for i, xyz in enumerate(xyz_vec):
                _curr_id += 1
                mons_list[i] = Atom("O", m=_m_mons)
                self.ionization_count += 1
                if (self.ionization_count % self.ionization_every == 0):
                    mons_list[i].q = -1*(abs(self.valency))
                mons_list[i].position = xyz
                mons_list[i].idx = _curr_id
                mons_list[i].ityp = 2
                mons_list[i].chainnum = chaincount
                if i < (npolychain - 2):
                    self.bonds.append(Bond(_curr_id, _curr_id + 1, bond_length))
            if (self.ionization_count % self.ionization_every == 0):
                mons_list[i].q = -1*(abs(self.valency))
            self.ionization_count += 1

            self.bonds.append(Bond(a1.idx, mons_list[0].idx, bond_length))
            self.bonds.append(Bond(a2.idx, mons_list[-1].idx, bond_length))
            self.mons.extend(mons_list)
            self.nmons = len(list(self.mons))
        self.number_of_chains = chaincount
        self.name_species.append('mons')
        self.species.append(self.mons)
        self.nspecies.append(self.nmons)
        self.mspecies.append(float(_m_mons))
        return None


    def create_cions(self):
        """
        generate the cions to have the total charge 0
        make sure you run this after you generate the gel
        # we will choose 1/6th of the atoms from the right side
        pick random ones and move them to the left
        """
        _m_cions = 1.
        shift_idx = self.nparticles
        print "shift_IDX for CIONS", shift_idx
        self.ncions = int(abs(self.charge)) # because now we have a valency
        self.cions = [None] * self.ncions
        R = np.random.uniform(0, 1, (self.ncions, 3))  * self.system_box
        for i, x_i in enumerate(R):
            self.cions[i] = Atom("H", m=_m_cions)
            self.cions[i].q = 1.
            self.cions[i].position = x_i
            self.cions[i].idx = shift_idx + i + 1
            self.cions[i].ityp = 3
        self.species.append(self.cions)
        self.nspecies.append(self.ncions)
        self.name_species.append("cions")
        self.mspecies.append(float(_m_cions))
        return None


    def create_probe_atoms(self):
        """
        generate a group of probe atoms
        that will be put in the system to probe the electrical charge
        """
        _m_probe = 1.
        shift_idx = self.nparticles
        print "shift_IDX for probe", shift_idx
        self.probes = [None] * self.nprobes
        R = np.random.uniform(0, 1, (self.nprobes, 3))  * self.system_box
        for i, x_i in enumerate(R):
            self.probes[i] = Atom("O", m=_m_probe)
            self.probes[i].q = 0.
            self.probes[i].position = x_i
            self.probes[i].idx = shift_idx + i + 1
            self.probes[i].ityp = 4
        self.species.append(self.probes)
        self.nspecies.append(self.nprobes)
        self.name_species.append("probes")
        self.mspecies.append(float(_m_probe))
        return None

    def create_dummy_probe(self):
        """creates a dummy probe atom so i can use them for my simulations"""
        _m_probe = 1.
        self.nspecies.append(self.nprobes)
        self.name_species.append("probes")
        self.mspecies.append(float(_m_probe))
        return None
    @staticmethod
    def n_from_c_salt(nnodes, nmons, ioniz, csalt):
        """gets number particles from the number of ionized gel ions, multiplied by the factor csalt
        nanions = (nnodes+nmons)*ioniz*csaltfactor
        """
        total_n_monomers = nnodes + nmons
        nsalt = int(csalt*total_n_monomers/ioniz)
        print "nsalt = ", nsalt
        return nsalt

    def create_salt(self):
        """
        generate the salt macro ions
        """
        self.nsalt = Mol.n_from_c_salt(self.nnodes, self.nmons, ioniz=self.ionization_every, csalt=self.csalt)
        # self.nsalt = Mol.n_from_c_salt(self.csalt, self.box[0])
        _m_cations = 1.0
        _m_anions = 1.0
        if self.nsalt > 0:
            # shift_idx = self.ncions + self.nnodes + self.nmons
            shift_idx = self.nparticles
            self.cations = [None] * (self.nsalt)
            self.anions = [None] * (self.nsalt)
            R = np.random.uniform(0, 1, (self.nsalt, 3)) * self.box
            for i, x_i in enumerate(R):
                self.cations[i] = Atom("Na",m=_m_cations)
                self.cations[i].q = 1.
                self.cations[i].position = x_i
                self.cations[i].idx = shift_idx + i + 1
                self.cations[i].ityp = 5
            shift_idx += self.nsalt
            # R = np.random.uniform(0, 1, (self.nsalt, 3)) * self.box
            R = np.random.rand(self.nsalt, 3) * self.box
            for i, x_i in enumerate(R):
                self.anions[i] = Atom("Cl",m=_m_anions)
                self.anions[i].q = -1.
                self.anions[i].position = x_i
                self.anions[i].idx = shift_idx + i + 1
                self.anions[i].ityp = 6
        self.name_species.append('cations')
        self.species.append(self.cations)
        self.nspecies.append(self.nsalt)
        self.mspecies.append(float(_m_cations))
        self.name_species.append('anions')
        self.species.append(self.anions)
        self.nspecies.append(self.nsalt)
        self.mspecies.append(float(_m_anions))
        return None



def write_data(mol, filename):
    """writes lammps data file"""
    outfile = (filename).rsplit('.', 1)[0] + '.data'
    nparticles = mol.nparticles
    nbond = mol.nbonds
    # xlow, xhigh = mol.cell_dims
    if mol.dummy_probe:
        ntypes = len(mol.name_species)
    else:
        ntypes = len(mol.species)
# -slab_thickness_dim*slab_l

    xlo = 0
    ylo = 0.
    zlo = 0.
    xhi = mol.system_box[0]
    zhi = mol.system_box[2]
    yhi = mol.system_box[1]

    # print "with types " + " ".join(str(list(mol.species)))
    print "with names " + " ".join(str(mol.name_species))
    print "with masses " + " ".join(str(mol.mspecies))

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
""".format(natoms=nparticles, nbonds=nbond,ntypes=ntypes, xlo=xlo, ylo=ylo, zlo=zlo, xhi=xhi, yhi=yhi, zhi=zhi)
# """.format(natoms=nparticles, nbonds=nbond,ntypes=ntypes, xlo=xlo, ylo=0., zlo=0., xhi=xhi, yhi=mol.system_box[1], zhi=mol.system_box[2])

    # print "mol.mspecies", mol.mspecies
    # print "mol.name_species", mol.name_species
    # print zip(mol.mspecies, mol.name_species)
    mass_string = ""
    for i, (ispec_mass, ispec_name) in enumerate(zip(mol.mspecies, mol.name_species)):
        mass_string += """{species_number} {species_mass} # {species_name}\n""".format(species_number=i+1, species_mass=ispec_mass, species_name=ispec_name)
    mass_string  += """\nAtoms\n"""

    with file(outfile, 'w') as f:
        f.write(header_string + "\n")
        f.write(mass_string + "\n")
        for atom_type, atom_type_name in zip(mol.species, mol.name_species):
            # atom_type_name = "lala"
            # print "type",atom_type
            # print "name",atom_type_name
            for ia, a in enumerate(atom_type):
                f.write("{index} {nmol} {ityp} {charge} {x} {y} {z} # {atom_type_name}\n".format(
                        index=a.idx, nmol=a.chainnum + 1, ityp=int(a.ityp),
                        charge=a.q, x=a.x, y=a.y, z=a.z, atom_type_name=atom_type_name))
        f.write('\nBonds\n\n')
        for ib, b in enumerate(mol.bonds):
            ib += 1
            # print b
            # f.write(" {0:7d} {1:7d}".format(b.j + 1, b.i + 1))
            f.write("{bond_id} {bond_type} {j:7d} {i:7d} \n".format(bond_id=ib,bond_type=b.ityp, j=b.j, i=b.i))
    return None

def read_cmd_line():
    """parses the command line parameters and return the args"""
    parser = argparse.ArgumentParser(description='program to create nodes, counter ions and salt ions (both cations and anions)')
    parser.add_argument('-f', '--filename', default='mol')
    parser.add_argument('-lt', '--lattice', default='diamond', help='generate the lattice of this type')
    parser.add_argument('-nm', '--nmonomers', default=40, type=int)
    parser.add_argument('-val', '--valency', default=1, type=int,help='valency of the backbone atoms')
    parser.add_argument('-l', '--length_scale', default=100., type=float, help='for diamond lattice the neighbor is ~0.4*scale, so for bondlength 1.0 one needs 4 monomers for length_scale=10.')
    parser.add_argument('--auto_length', dest='auto_length', action='store_true', help='automatically calculate the chain length based on the equilibrium length of the Flory chain in a good solvent Rnu = bN^nu, nu = 0.588, for diamond lattice the neighbor is Rnu = l = N^0.588/0.4')
    parser.add_argument('--rw', dest='random_walk_chain', action='store_true', help='connect the nodes via a von mises random walk - good for dense systems, otherwise a straight chain is used, which needs to be relaxed separately. use with care since it doesnt always work well. use straight chain when using no salt')
    parser.add_argument('--switch_xz', dest='switch_xz', action='store_true', help='switch between x and z coordinates, so lammps can use kspace slab')
    parser.add_argument('-a', '--alpha', default=1., type=int, help='degree of ionization, only 1,alpha+1,2alpha+1.... will be charged, 1 - half of the monomers will be charged')
    parser.add_argument('-df', '--density_factor', default=1., type=float, help='degree of difference, by giving the chain length scaling factor, chains on the left are nm*densityfactor long on the right nm long')
    # parser.add_argument('-s', '--salt', default=0., type=float, help='concentration of salt ions, no macro-ions if 0.')
    parser.add_argument('-s', '--salt', default=0, type=float, help='do you want to have salt ions in your system?  If 0 - then no salt, if 1 then amount of salt anions will be equal to the amount of ionized gel ions in the system, 2 - then twice')
    # parser.add_argument('-d', '--dimensions', default=[3, 3, 3])
    parser.add_argument('-d','--dimensions', nargs='+', type=int, default=[1, 1,1])
    parser.add_argument('-p','--pbc', nargs='+', type=int, default=[1, 1, 1])
    args = parser.parse_args()
    if args.auto_length == True:
        print "doing autolength"
        scale = float(args.density_factor)
        n = int(args.nmonomers)
        nscaled = n * scale
        l_flory = 1. * nscaled**0.588
        print "flory length = ", l_flory
        if args.lattice == 'sc':
            args.length_scale=l_flory
            print "length set to flory_length (for sc neigh distnace) = ", l_flory
        else:
            args.length_scale=l_flory/.4
            print "length set to flory_length/.4 (for diamond neigh distnace) = ", l_flory/.4
    else:
         print "nothing done automatically"
    print"Check that args.do_something=" + str(args.auto_length) + " is always a bool"
    return args

#vmd
def write_vmd(filename, length_scale, mol):
    """writes the vmd file for vizualization"""

    input_data = (filename).rsplit('.', 1)[0] + '.data'
    viz_file = (filename).rsplit('.', 1)[0] + '.vmd'
    size_dict = dict()
    size_dict['mons'] = length_scale / 36. # r_cpk
    size_dict['nodes'] = length_scale / 5. # r_vdw
    size_dict['anions'] = length_scale / 30. # r_salt
    size_dict['cations'] = length_scale / 30. # r_salt
    size_dict['cions'] = 3.5  # r_point
    size_dict['probes'] = 2. # r_probe

    viz_dict = dict()
    viz_dict['mons'] = "VDW" # r_cpk
    viz_dict['nodes'] = "VDW" # r_vdw
    viz_dict['anions'] = "VDW" # r_salt
    viz_dict['cations'] = "VDW" # r_salt
    viz_dict['cions'] = "Points"  # r_point
    viz_dict['probes'] = "VDW" # r_probe

    col_dict = dict()
    col_dict['mons'] = 0 # r_cpk
    col_dict['nodes'] = 1 # r_vdw
    col_dict['anions'] = 27 # r_salt
    col_dict['cations'] = 16 # r_salt
    col_dict['cions'] = 7  # r_point
    col_dict['probes'] = 14 # r_probe

    header_string="""\
color Display Background white
topo readlammpsdata {input_data} full waitfor all""".format(input_data=input_data)

    middle_string = ""
    for i, (atom_type, atom_type_name) in enumerate(zip(mol.species, mol.name_species)):
        if i > 0:
            middle_string += """mol addrep 0"""
        type_idx = atom_type[0].ityp
        if viz_dict[atom_type_name] == 'Points':
            resolution = ""
        else:
            resolution = "12.0"
        middle_string += """\n
            # {atom_type_name}
            mol modselect {ii} 0 name is "{ityp}"
            mol modcolor {ii} 0 ColorID {col_type}
            mol modstyle {ii} 0 {viz_type} {r} {resolution} \n
            """.format(ii=i, ityp=atom_type_name, col_type=col_dict[atom_type_name],viz_type=viz_dict[atom_type_name],r=size_dict[atom_type_name], atom_type_name=atom_type_name, resolution=resolution)
    footer_string = """
        set sel_all [atomselect top "name is cions"]
        set minmax [measure minmax $sel_all]
        set box [vecscale 1.0 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
        pbc set $box
        pbc box
        display resetview
        # render TachyonInternal {input_data}.png /usr/bin/open %s
    """.format(input_data=input_data)

    with file(viz_file, 'w') as f:
        f.write(header_string+middle_string+footer_string)

def main():
    """program to generate the  pe lattice that can be used as a module"""
    args = read_cmd_line()
    verbose = False
    print args
    Gel = Mol(topol=args.lattice, scale=args.length_scale, npoly=args.nmonomers, n_=args.dimensions, ionization_every=args.alpha, csalt=args.salt,densityfactor=args.density_factor, valency=args.valency, pbc_=args.pbc, rand_walk_chain=args.random_walk_chain)
    if verbose:
        print Gel
    print Gel
    print Gel.nparticles
    # print Gel.bottom_sites
    # print Gel.top_sites
    write_data(Gel, args.filename)
    # write_vmd(args.filename, args.length_scale)
    write_vmd(args.filename, 10, Gel)

    # print "System, nodes = {0}, mons = {1}, neighs = {2}".format(Gel.nnodes, Gel.nmons, Gel.nneighs)
if __name__ == '__main__':
    main()
