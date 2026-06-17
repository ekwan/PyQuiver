"""Parser for Gaussian output files (g09 / g16)."""

import os
import re
import logging

import numpy as np

from ._common import ParsedSystem, parse_serial_lower_hessian

logger = logging.getLogger("pyquiver")


def _tail(filename):
    """Return the last line of a file (used to check normal termination)."""
    # https://stackoverflow.com/questions/46258499
    with open(filename, 'rb') as f:
        try:  # catch OSError in case of a one-line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode()


def _valid_geom_line(split_line):
    if len(split_line) == 6:
        try:
            int(split_line[0]); int(split_line[1]); int(split_line[2])
            return True
        except ValueError:
            return False
    return False


def _parse_hessian(data, filename):
    # the force constant matrix lives in the archive between 1\1\GINC and \@
    raw_archive = re.findall(r"1\\1\\GINC(.+?)\\(\s*)@", data, re.DOTALL)
    archive = None
    for candidate in raw_archive:
        squashed = re.sub(r'[\s+]', '', candidate[0])
        archive = re.search(r"NImag\=(.+?)$", squashed, re.DOTALL)
        if archive:
            break
    if not archive:
        raise ValueError("No frequency job detected in %s." % filename)

    return archive.group(0).split('\\')[2].split(',')


def parse(path):
    """Parse a Gaussian output file into a :class:`ParsedSystem`."""
    # require a normally terminated job (snippet files are exempt)
    if not path.endswith(".snip") and "Normal termination" not in _tail(path):
        raise ValueError("Gaussian job %s terminated in an error" % path)

    with open(path, 'r') as f:
        data = f.read()

    # the verbose route card (#p) is required for the archive to be present
    if re.search(r" *#[pP] ", data) is None and not path.lower().endswith(".snip"):
        raise ValueError(
            "Gaussian output file %s was not run with the verbose flag, so it "
            "does not contain enough information for PyQuiver to run. Please "
            "re-run this calculation with a route card that starts with #p"
            % path)

    m = re.search(r"NAtoms\= +([0-9]+)", data)
    if not m:
        raise ValueError("Number of atoms not detected.")
    number_of_atoms = int(m.group(1))

    atomic_numbers = [0] * number_of_atoms
    positions = np.zeros((number_of_atoms, 3))

    # prefer standard orientation; fall back to input orientation (nosymm)
    m = None
    for m in re.finditer(r"Standard orientation(.+?)Rotational constants \(GHZ\)",
                         data, re.DOTALL):
        pass
    if m is None:
        for m in re.finditer(r"Input orientation(.+?)Distance matrix",
                             data, re.DOTALL):
            pass
        if m is not None:
            logger.info("Couldn't find standard orientation so used input orientation instead.")
    if m is None:
        for m in re.finditer(r"Input orientation(.+?)Rotational constants \(GHZ\)",
                             data, re.DOTALL):
            pass
        if m is not None:
            logger.info("Couldn't find standard orientation so used input orientation instead.")
    if m is None:
        raise ValueError("Geometry table not detected.")

    for line in m.group(1).split('\n'):
        fields = [x for x in line.split(' ') if x]
        if _valid_geom_line(fields):
            center = int(fields[0]) - 1
            atomic_numbers[center] = int(fields[1])
            for e in range(3):
                positions[center][e] = fields[3 + e]

    raw_fcm = _parse_hessian(data, path)
    hessian = parse_serial_lower_hessian(raw_fcm, number_of_atoms)

    return ParsedSystem(atomic_numbers, positions, hessian)
