import os
import logging

import numpy as np

from .constants import PHYSICAL_CONSTANTS, LINEARITY_THRESHOLD, DROP_NUM_LINEAR
from . import parsers
from .parsers import native

logger = logging.getLogger("pyquiver")

# represents a geometric arrangement of atoms with specific masses
class Isotopologue(object):
    def __init__(self, id_, system, masses):
        self.name = id_
        self.system = system
        self.masses = masses
        self.frequencies = None

        self.number_of_atoms = system.number_of_atoms
        self.mw_hessian = self.calculate_mw_hessian()

    def __str__(self):
        returnString  = "Isotopologue: %s\n" % self.name
        returnString += " masses: %s" % self.masses.__str__()
        return returnString

    def calculate_mw_hessian(self):
        # mass-weight the Hessian: H'_ij = H_ij / sqrt(m_i m_j), where each
        # atom's mass is repeated across its three Cartesian coordinates
        masses_ij = np.outer(np.array(self.masses), np.array(self.masses))
        inv_sqrt_masses = 1 / np.sqrt(masses_ij)
        inv_sqrt_masses = np.repeat(inv_sqrt_masses, 3, 0)
        inv_sqrt_masses = np.repeat(inv_sqrt_masses, 3, 1)
        return self.system.hessian * inv_sqrt_masses

    def calculate_frequencies(self, imag_threshold, scaling=1.0, method="mass weighted hessian"):
        # short circuit if frequencies have already been calculated
        if self.frequencies is not None:
            return self.frequencies

        if method == "mass weighted hessian":
            imaginary_freqs = []
            small_freqs = []
            conv_factor = PHYSICAL_CONSTANTS['Eh']/(PHYSICAL_CONSTANTS['a0']**2 * PHYSICAL_CONSTANTS['amu'])

            v = np.linalg.eigvalsh(self.mw_hessian*conv_factor)

            constant = scaling / (2*np.pi*PHYSICAL_CONSTANTS['c'])
            freqs = np.sqrt(np.abs(v)) * np.sign(v) * constant

            #freqs2 = [ np.copysign(np.sqrt(np.abs(freq)),freq) * constant for freq in v ]
            #assert np.allclose(np.array(freqs), freqs2)

            freqs.sort()

            imaginary_freqs = []
            small_freqs = []
            regular_freqs = []

            # detect imaginary frequencies
            for f in freqs:
                if f < -imag_threshold:
                    imaginary_freqs.append(f)

            if len(imaginary_freqs) > 1:
                logger.warning("multiple imaginary frequencies detected")

            # strip the imaginary frequencies
            freqs = freqs[len(imaginary_freqs):]

            if self.system.is_linear:
                small_freqs = freqs[:DROP_NUM_LINEAR]
                regular_freqs = freqs[DROP_NUM_LINEAR:]
            else:
                small_freqs = freqs[:1+DROP_NUM_LINEAR]
                regular_freqs = freqs[1+DROP_NUM_LINEAR:]

            # bugfix 2/6/20: third argument is regular_freqs, not freqs!
            self.frequencies = (small_freqs, imaginary_freqs, np.array(regular_freqs), len(small_freqs))
            return self.frequencies
        else:
            raise ValueError("unknown frequency calculation type")


class System(object):
    # represents a molecule's geometry and Cartesian Hessian, read from an
    # electronic-structure output file via the parsers package
    def __init__(self, outfile, style="gaussian"):
        self.filename = outfile
        self.is_linear = True

        logger.info("Reading data from %s with style %s", outfile, style)

        parsed = parsers.parse(outfile, style)
        self.atomic_numbers = parsed.atomic_numbers
        self.number_of_atoms = len(parsed.atomic_numbers)
        self.hessian = parsed.hessian
        self.positions_angstrom = parsed.positions_angstrom
        # geometry in bohr, as needed for the mass-weighted analysis
        self.positions = parsed.positions_angstrom / PHYSICAL_CONSTANTS['atb']

        self._detect_linear()

    def _detect_linear(self):
        # take every pair of bonds sharing atom 0; if any pair is sufficiently
        # non-parallel the molecule is non-linear
        positions = self.positions_angstrom
        if self.number_of_atoms > 2:
            for i in range(1, len(positions) - 1):
                for j in range(i + 1, len(positions)):
                    diff0i = positions[0] - positions[i]
                    diff0j = positions[0] - positions[j]
                    u = diff0i / np.linalg.norm(diff0i)
                    v = diff0j / np.linalg.norm(diff0j)
                    if 1 - np.abs(np.inner(u, v)) > LINEARITY_THRESHOLD:
                        self.is_linear = False
                        break
                if not self.is_linear:
                    break
        logger.debug("Molecule is %s.", "linear" if self.is_linear else "not linear")

    def dump_pyquiver_input_file(self, extension=".qin"):
        path = os.path.splitext(self.filename)[0] + extension
        serial = native.serialize(self)
        with open(path, 'w') as f:
            f.write(serial)
        return serial
