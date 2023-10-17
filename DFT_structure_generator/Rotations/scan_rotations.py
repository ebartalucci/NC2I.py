import numpy as np
import os
import argparse
import math
from copy import deepcopy


class Atom:
    def __init__(self, num, frag, line):
        self.number = num
        self.fragment = frag
        self.symbol = line.split()[0]
        self.x = float(line.split()[1])
        self.y = float(line.split()[2])
        self.z = float(line.split()[3])

    def get_atom_line(self):
        return "{:2} {:12.7f} {:12.7f} {:12.7f}".format(self.symbol, self.x, self.y, self.z)

    def move_atom(self, disp):
        self.x += disp['x']
        self.y += disp['y']
        self.z += disp['z']


def get_distance(point1, point2):
    return math.sqrt((point1['x'] - point2['x'])**2 + (point1['y'] - point2['y'])**2 + (point1['z'] - point2['z'])**2)


def get_angle(v1, v2):
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    cos_angle = np.dot(v1, v2) / (v1_norm * v2_norm)
    angle = math.degrees(math.acos(cos_angle))
    return angle


def read_settings(path):
    with open(path, 'r') as inp:
        sets = inp.readlines()
        settings = {}
        for index, line in enumerate(sets):
            if '$fragment1' in line:
                settings['fragment1'] = [int(at) for at in sets[index+1].split()]
            if '$fragment2' in line:
                settings['fragment2'] = [int(at) for at in sets[index+1].split()]
            if '$fixpoint1' in line:
                tmp = sets[index+1].split()
                settings['fixpoint1'] = {'x': float(tmp[0]), 'y': float(tmp[1]), 'z': float(tmp[2])}
            if '$fixpoint2' in line:
                tmp = sets[index+1].split()
                settings['fixpoint2'] = {'x': float(tmp[0]), 'y': float(tmp[1]), 'z': float(tmp[2])}
            if '$angles' in line:
                tmp = sets[index+1].split()
                settings['angle_start'] = int(tmp[0])
                settings['angle_end'] = int(tmp[1])
                settings['angle_step'] = int(tmp[2])
    return settings


def read_atoms(path, sets):
    with open(path, 'r') as inp:
        data = inp.readlines()
        nat = int(data[0].split()[0])
        atoms = []
        for index, line in enumerate(data[2:], start=2):
            atnum = index - 1
            if atnum in sets['fragment1']:
                atfrag = 1
            elif atnum in sets['fragment2']:
                atfrag = 2
            else:
                print("ERROR: Atom {} has not been assigned to fragment 1 or 2!".format(atnum))
                exit()
            atoms.append(Atom(atnum, atfrag, line))

    if nat != len(atoms):
        print("ERROR: Number of atoms is not equal to the elements in list atoms!")
        exit()

    return atoms


def get_molecule_lines(atom_list, comment=""):
    mol = []
    mol.append(str(len(atom_list)))
    mol.append(comment)
    for at in atom_list:
        mol.append(at.get_atom_line())
    return mol


def print_molecule(atom_list, comment=""):
    mol = get_molecule_lines(atom_list, comment)
    for line in mol:
        print(line)


def write_molecule(path, outname, atom_list, comment=""):
    mol = get_molecule_lines(atom_list, comment)
    with open(os.path.join(path, outname), 'w') as out:
        out.write("\n".join(mol) + "\n")


def get_new_fixpoint(fix1, fix2, angle):
    v1 = np.array([fix1['x'], fix1['y'], fix1['z']])
    v2 = np.array([fix2['x'], fix2['y'], fix2['z']])
    axis = np.cross(v1, v2)
    axis /= np.linalg.norm(axis)
    rot_mat = rotation_matrix(axis, math.radians(angle))
    v2_rotated = np.dot(rot_mat, v2)
    fix_new = {
        'x': v2_rotated[0],
        'y': v2_rotated[1],
        'z': v2_rotated[2]
    }
    return fix_new


def rotation_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])


parser = argparse.ArgumentParser()
parser.add_argument('structure', type=str, help='input structure file in .xyz format', metavar='struct')
parser.add_argument('settings', type=str, help='input settings in the given format in some file like scan.inp', metavar='sets')
args = parser.parse_args()

workdir = os.getcwd()
inputstructure = args.structure
inputsettings = args.settings

settings = read_settings(os.path.join(workdir, inputsettings))
atoms = read_atoms(os.path.join(workdir, inputstructure), settings)

fixpoint_angle = get_angle(
    [settings['fixpoint1']['x'], settings['fixpoint1']['y'], settings['fixpoint1']['z']],
    [settings['fixpoint2']['x'], settings['fixpoint2']['y'], settings['fixpoint2']['z']]
)

desired_angles = range(settings['angle_start'], settings['angle_end'], settings['angle_step'])
desired_angles.append(settings['angle_end'])

for index, angle in enumerate(desired_angles):
    name = inputstructure[:-4]
    suffix = "_" + str(index+1).zfill(3) + "_" + str(angle) + ".xyz"

    fixpoint_new = get_new_fixpoint(settings['fixpoint1'], settings['fixpoint2'], angle)
    disp = {
        'x': fixpoint_new['x'] - settings['fixpoint2']['x'],
        'y': fixpoint_new['y'] - settings['fixpoint2']['y'],
        'z': fixpoint_new['z'] - settings['fixpoint2']['z']
    }

    atoms_moved = deepcopy(atoms)
    for atom in atoms_moved:
        atom.move_atom(disp)

    write_molecule(workdir, name + suffix, atoms_moved)
