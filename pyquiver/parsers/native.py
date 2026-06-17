"""Parser and serializer for PyQuiver's native ``.qin`` format.

The format is::

    <number of atoms>
    <center>,<atomic number>,<x>,<y>,<z>      (one line per atom, angstroms)
    <comma-separated lower-triangular Hessian>
"""

import numpy as np

from ._common import ParsedSystem, parse_serial_lower_hessian


def parse(path):
    """Parse a ``.qin`` file into a :class:`ParsedSystem`."""
    with open(path, 'r') as f:
        lines = f.read().split("\n")

    try:
        number_of_atoms = int(lines[0])
    except ValueError:
        raise ValueError("first line must contain integer number of atoms.")

    atomic_numbers = [0] * number_of_atoms
    positions = np.zeros((number_of_atoms, 3))
    for line in lines[1:number_of_atoms + 1]:
        fields = line.split(',')
        try:
            center_number, atomic_number, x, y, z = fields
        except ValueError:
            raise ValueError("the following line in the geometry did not have "
                             "the appropriate number of fields: {0}".format(line))
        center = int(center_number)
        atomic_numbers[center] = int(atomic_number)
        positions[center][0] = x
        positions[center][1] = y
        positions[center][2] = z

    fcm_fields = lines[number_of_atoms + 1].split(',')
    hessian = parse_serial_lower_hessian(fcm_fields, number_of_atoms)

    return ParsedSystem(atomic_numbers, positions, hessian)


def serialize(system):
    """Serialize a System to native ``.qin`` text."""
    serial = "%d\n" % system.number_of_atoms
    for i in range(system.number_of_atoms):
        serial += "{0},{1},{2},{3},{4}\n".format(
            i, system.atomic_numbers[i],
            system.positions_angstrom[i, 0],
            system.positions_angstrom[i, 1],
            system.positions_angstrom[i, 2])
    for i in range(system.number_of_atoms * 3):
        for j in range(0, i + 1):
            serial += str(system.hessian[i, j]) + ","
    return serial
