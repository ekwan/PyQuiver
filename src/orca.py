import re
import numpy as np

from io import StringIO

DEBUG = False

elemToNum = {
    "H": 1, "He": 2,
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38,
    "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56,
    "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88,
    "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112
}


def debug(s):
    if DEBUG:
        print("DEBUG: {}".format(s))

def bohr2angstrom(length):
    return length*0.529177249

def parse_orca_output(out_data):
    lines = out_data.split('\n')

    # Molecule position data
    try:
        atoms_line = next(
            (n for n, line in enumerate(lines)
             if line.startswith('$atoms')))
    except StopIteration:
        raise ValueError("Could not find '$atoms' in output data.")

    number_of_atoms = int(lines[atoms_line+1])
    debug("number_of_atoms={}".format(number_of_atoms))

    atomic_numbers = []
    positions = np.zeros(shape=(number_of_atoms, 3))
    for i, line in enumerate(
            lines[atoms_line+2:atoms_line+number_of_atoms+2]):
        _, atom, _, x, y, z = re.split(r'\s+', line)

        atomic_numbers.append(elemToNum[atom])
        positions[i][0] = bohr2angstrom(float(x))
        positions[i][1] = bohr2angstrom(float(y))
        positions[i][2] = bohr2angstrom(float(z))

    debug("atomic_numbers={}".format(atomic_numbers))
    debug("positions={}".format(positions))

    try:
        hessian_line = next(
            (n for n, line in enumerate(lines)
             if line.startswith('$hessian')))
    except StopIteration:
        raise ValueError("Could not find '$hessian' in output data.")

    hessian_size = int(lines[hessian_line+1])
    debug("Size of hessian: {}".format(hessian_size))

    hessian = None

    for i, batch in enumerate(range(0, hessian_size, 5)):
        batch_size = min(5, hessian_size-batch)

        s = StringIO(u"\n".join(
            lines[hessian_line+(i*(hessian_size+1)+2):
                  hessian_line+((i+1)*(hessian_size+1))+2]))

        h = np.loadtxt(s, skiprows=1,
                       converters=dict((i+1, float) for i in range(batch_size)),
                       usecols=tuple(i+1 for i in range(batch_size)))

        hessian = h if hessian is None else np.hstack((hessian, h))

    hessian = (hessian + hessian.T)/2
    debug("Hessian: {}".format(hessian))

    return atomic_numbers, positions, hessian


if __name__ == '__main__':
    DEBUG = True
    import sys
    with open(sys.argv[1], 'r') as f:
        out_data = f.read()
    parse_orca_output(out_data)
