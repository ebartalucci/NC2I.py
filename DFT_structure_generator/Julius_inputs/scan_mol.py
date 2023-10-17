#!/usr/bin/env python3


import os
import argparse
import math
from copy import deepcopy


# atom class stores all essential information about each atom
class atom:
    def __init__(self, num, frag, line):
        self.number = num                 # atom number in input .xyz file
        self.fragment = frag              # atom belongs to fragment frag (1 or 2)
        self.symbol = line.split()[0]     # element symbol
        self.x = float(line.split()[1])   # x-coordinate (angstrom)
        self.y = float(line.split()[2])   # y-coordinate (angstrom)
        self.z = float(line.split()[3])   # z-coordinate (angstrom)

    # get string line for this atom in .xyz file format
    def get_atom_line(self):
        return "{:2} {:12.7f} {:12.7f} {:12.7f}".format(self.symbol, self.x, self.y, self.z)

    # move the atom by the displacement vector disp by adding x,y,z coordinates
    def move_atom(self, disp):
        self.x += disp['x']
        self.y += disp['y']
        self.z += disp['z']


# get distance (angstrom) between two points in dicts
def get_distance(point1, point2):
    return math.sqrt((point1['x']-point2['x'])**2 + (point1['y']-point2['y'])**2 + (point1['z']-point2['z'])**2)


# read all data from the input settings file and store in a dict
def read_settings(path):
    with open(path, 'r') as inp:
        sets = inp.readlines()
        settings = {}
        for index, line in enumerate(sets):
            if '$fragment1' in line:
                settings['fragment1'] = [int(at) for at in sets[index+1].split()]   # list of atoms that belong to fragment 1
            if '$fragment2' in line:
                settings['fragment2'] = [int(at) for at in sets[index+1].split()]   # list of atoms that belong to fragment 2
            if '$fixpoint1' in line:
                tmp = sets[index+1].split()
                settings['fixpoint1'] = {'x': float(tmp[0]), 'y': float(tmp[1]), 'z': float(tmp[2])}   # fixpoint of fragment 1
            if '$fixpoint2' in line:
                tmp = sets[index+1].split()
                settings['fixpoint2'] = {'x': float(tmp[0]), 'y': float(tmp[1]), 'z': float(tmp[2])}   # fixpoint of fragment 2
            if '$distances' in line:
                tmp = sets[index+1].split()
                settings['dist_start'] = int(tmp[0])
                settings['dist_end'] = int(tmp[1])
                settings['dist_step'] = int(tmp[2])
    return settings


# read the molecule .xyz file and store information in a list of atom objects
def read_atoms(path, sets):
    with open(path, 'r') as inp:
        data = inp.readlines()
        nat = int(data[0].split()[0])                 # number of atoms
        atoms = []
        for index, line in enumerate(data[2:], start=2):
            atnum = index - 1                         # index starts at 0, but atoms start in 3rd line
            if atnum in sets['fragment1']: atfrag = 1
            elif atnum in sets['fragment2']: atfrag = 2
            else:
                print("ERROR: Atom {} has not been assigned to fragment 1 or 2!".format(atnum))
                exit()
            atoms.append(atom(atnum, atfrag, line))   # store atom objects in list atoms
    
    # check if nat is the same as the number of elements in list atoms
    if nat != len(atoms):
        print("ERROR: Number of atoms is not equal to the elements in list atoms!")
        exit()

    return atoms


# get the whole molecule as list of string lines in .xyz format
# requires a list atoms of objects atom
# can print an optional comment in the second line
def get_molecule_lines(atom_list, comment=""):
    mol = []
    mol.append(str(len(atom_list)))
    mol.append(comment)
    for at in atom_list:
        mol.append(at.get_atom_line())
    return mol


# print the molecule in .xyz format from get_molecule_lines
def print_molecule(atom_list, comment=""):
    mol = get_molecule_lines(atom_list, comment)
    for line in mol: print(line)


# write the molecule in .xyz file in path/outname
def write_molecule(path, outname, atom_list, comment=""):
    mol = get_molecule_lines(atom_list, comment)
    with open(os.path.join(path, outname), 'w') as out:
        out.write("\n".join(mol) + "\n")


# get a point 2' on the line l(d) spanned by fixpoints 1 and 2 in dependence on distance d
# d: distance of the new point 2' to fixpoint 2 behind fixpoint 1: --1----2--2'--
# r_1, r_2, r_2': vectors of the respective points 1, 2, and 2'
# l(d) = r_2' = r_2 + (r_2-r_1)/|r_2-r_1| * d
# fix1: r_1, fix2: r_2, dist: d (desired distance), d_12: distance between r_1 and r_2, fix_new: r_2'
def get_new_fixpoint(fix1, fix2, dist):
    d_12 = get_distance(fix1, fix2)
    fix_new = {
        'x': fix2['x'] + ((fix2['x']-fix1['x'])/d_12) * dist,
        'y': fix2['y'] + ((fix2['y']-fix1['y'])/d_12) * dist,
        'z': fix2['z'] + ((fix2['z']-fix1['z'])/d_12) * dist
    }
    return fix_new


# from an old and a new fixpoint, get displacement vector that maps old on new: r_new - r_old
def get_displacement_vector(fix_old, fix_new):
    return {'x': fix_new['x']-fix_old['x'], 'y': fix_new['y']-fix_old['y'], 'z': fix_new['z']-fix_old['z']}


# set up parser and parse the input files
parser = argparse.ArgumentParser()
parser.add_argument('structure', type=str, help='input structure file in .xyz format', metavar='struct')
parser.add_argument('settings', type=str, help='input settings in the given format in some file like scan.inp', metavar='sets')
args = parser.parse_args()

# get environment
workdir = os.getcwd()
inputstructure = args.structure
inputsettings = args.settings

# read the input settings and geometry from the files
settings = read_settings(os.path.join(workdir, inputsettings))
atoms = read_atoms(os.path.join(workdir, inputstructure), settings)

# print input molecule (for testing reasons)
#print_molecule(atoms, "old")

# get the original fixpoint distance
fixpoint_distance = get_distance(settings['fixpoint1'], settings['fixpoint2'])

# get the distances relative to the fixpoint distance from the scan.inp file
# e.g. -50 200 10: -50% to 200% in 10% steps of the fixpoint distance
desired_distances = [(i/100.0) for i in range(settings['dist_start'], settings['dist_end'], settings['dist_step'])]
desired_distances.append(settings['dist_end']/100.0)

# generate new files with the desired distances
for index, dist in enumerate(desired_distances):
    name = inputstructure[:-4]                                                 # removes last 4 characters from input file name (should be ".xyz")
    suffix = "_" + str(index+1).zfill(3) + "_" + str(int(dist*100)) + ".xyz"   # adds _number_step%.xyz zu the file name

    # get new fixpoint and displacement vector
    fixpoint_new = get_new_fixpoint(settings['fixpoint1'], settings['fixpoint2'], dist*fixpoint_distance)
    disp = get_displacement_vector(settings['fixpoint2'], fixpoint_new)

    # get new molecule where fragment 2 is moved according to the displacement vector
    atoms_moved = deepcopy(atoms)
    for at in atoms_moved:
        if at.fragment == 2: at.move_atom(disp)

    # write the new .xyz file
    write_molecule(workdir, name + suffix, atoms_moved, inputstructure + " with distance " + str(int(dist*100)) + "%")
