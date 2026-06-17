"""Parser for ORCA .hess files."""

import re
import logging
from io import StringIO

import numpy as np

from ..constants import SYMBOL_TO_Z, PHYSICAL_CONSTANTS
from ._common import ParsedSystem

logger = logging.getLogger("pyquiver")

_BOHR_TO_ANGSTROM = PHYSICAL_CONSTANTS["atb"]


def parse_orca_output(out_data):
    """Parse the text of an ORCA .hess file.

    Returns ``(atomic_numbers, positions_angstrom, hessian)``.
    """
    lines = out_data.split('\n')

    # geometry
    try:
        atoms_line = next(n for n, line in enumerate(lines)
                          if line.startswith('$atoms'))
    except StopIteration:
        raise ValueError("Could not find '$atoms' in output data.")

    number_of_atoms = int(lines[atoms_line + 1])

    atomic_numbers = []
    positions = np.zeros((number_of_atoms, 3))
    for i, line in enumerate(lines[atoms_line + 2:atoms_line + number_of_atoms + 2]):
        _, atom, _, x, y, z = re.split(r'\s+', line)
        atomic_numbers.append(SYMBOL_TO_Z[atom])
        positions[i][0] = float(x) * _BOHR_TO_ANGSTROM
        positions[i][1] = float(y) * _BOHR_TO_ANGSTROM
        positions[i][2] = float(z) * _BOHR_TO_ANGSTROM

    # hessian
    try:
        hessian_line = next(n for n, line in enumerate(lines)
                            if line.startswith('$hessian'))
    except StopIteration:
        raise ValueError("Could not find '$hessian' in output data.")

    hessian_size = int(lines[hessian_line + 1])
    hessian = None
    for i, batch in enumerate(range(0, hessian_size, 5)):
        batch_size = min(5, hessian_size - batch)
        s = StringIO(u"\n".join(
            lines[hessian_line + (i * (hessian_size + 1) + 2):
                  hessian_line + ((i + 1) * (hessian_size + 1)) + 2]))
        h = np.loadtxt(s, skiprows=1,
                       converters=dict((j + 1, float) for j in range(batch_size)),
                       usecols=tuple(j + 1 for j in range(batch_size)),
                       ndmin=2)
        hessian = h if hessian is None else np.hstack((hessian, h))

    hessian = (hessian + hessian.T) / 2
    return atomic_numbers, positions, hessian


def parse(path):
    """Parse an ORCA .hess file at ``path`` into a :class:`ParsedSystem`."""
    with open(path, 'r') as f:
        atomic_numbers, positions, hessian = parse_orca_output(f.read())
    return ParsedSystem(atomic_numbers, positions, hessian)
