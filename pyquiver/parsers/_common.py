"""Shared helpers and the common return type for input-file parsers."""

from collections import namedtuple

import numpy as np

# the normalized result every parser returns:
#   atomic_numbers     - list[int], one per atom
#   positions_angstrom - np.ndarray, shape (n_atoms, 3)
#   hessian            - np.ndarray, shape (3*n_atoms, 3*n_atoms), in hartree/bohr^2
ParsedSystem = namedtuple("ParsedSystem",
                          ["atomic_numbers", "positions_angstrom", "hessian"])


def parse_serial_lower_hessian(fields, number_of_atoms):
    """Expand a serialized lower-triangular Hessian into the full matrix.

    ``fields`` is the comma-split lower triangle (row-major). A trailing comma
    in the source leaves an empty final field, so blanks are dropped before
    converting to float.
    """
    fields = [f for f in fields if f.strip() != ""]

    # map every (i, j) of the full grid onto its lower-triangle serial index
    range1d = np.arange(3 * number_of_atoms, dtype=int)
    xgrid, ygrid = np.meshgrid(range1d, range1d)
    agrid = np.minimum(xgrid, ygrid)
    bgrid = np.maximum(xgrid, ygrid)
    idxs = (bgrid * (bgrid + 1) / 2 + agrid).astype(int)

    return np.array(fields, dtype=float)[idxs]
