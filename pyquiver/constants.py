# This file holds physical constants and reads atomic weights.
import os
import re
import logging

logger = logging.getLogger("pyquiver")
###############

# Physical Constants

PHYSICAL_CONSTANTS = {
    'h'  : 6.626070E-34, # Planck's constants in J * s
    'c'  : 2.997925E+10, # speed of light in units of cm/s
    'Eh' : 4.359745E-18, # energy of a hartree in units of J = kg m^2/s^2
    'a0' : 5.291772E-11, # bohr radius in m
    'atb': 5.291772E-01, # angstroms per bohr
    'amu': 1.660468E-27, # atomic mass unit in units kg
    'kB' : 1.380649E-23  # Boltzmann's constant in J/K
}
    #CM/2.998E10/,EM/1.440E13/,HBC/1.4387/

###############

# Atomic Weight Information

class Element(object):
    def __init__(self, full_name, atomic_number, symbol, default_mass):
        # the name of this element, like "hydrogen"
        full_name = str(full_name)
        self.full_name = full_name
        if re.match("[^a-z]", full_name):
            raise ValueError("Unexpected non-lowercase character in element "
                             "name: %s" % full_name)

        # the symbol of this element, like "H"
        symbol = str(symbol)
        self.symbol = symbol
        if re.match("[^a-zA-Z]", symbol):
            raise ValueError("Unexpected non-letter character in element "
                             "symbol: %s" % symbol)
        if len(symbol) < 1 or len(symbol) > 2:
            raise ValueError("Unexpected length of element symbol (must be 1 "
                             "or 2): %s" % symbol)

        # the atomic number of this element, like 1
        atomic_number = int(atomic_number)
        self.atomic_number = atomic_number
        if atomic_number < 1 or atomic_number > 200:
            raise ValueError("Unexpected atomic number: %d" % atomic_number)

        # the average weight for this element, like 1.00783
        default_mass = float(default_mass)
        self.default_mass = default_mass
        if default_mass < 0.0 or default_mass > 500.0:
            raise ValueError("Unexpected default mass: %f" % default_mass)

        # pairs of tuples strings (like "2H") to masses (like 2.0141)
        self.replacements = []

    def __str__(self):
        string = "%s (%s, Z=%d, default mass = %.4f" % (self.full_name.capitalize(), self.symbol, self.atomic_number, self.default_mass)
        if len(self.replacements) == 0:
            string += ", no isotopic replacements possible)\n"
        else:
            string += ")\n"
            for s,m in self.replacements:
                string += "    %2s : %.4f\n" % (s,m)
        return string[:-1]

    def add_replacement(self, symbol, mass):
        symbol = str(symbol)
        if re.match("[^a-zA-Z0-9]", symbol):
            raise ValueError("Unexpected non-letter character in isotopic "
                             "replacement symbol: %s" % symbol)
        if len(symbol) < 1 or len(symbol) > 4:
            raise ValueError("Unexpected length of element symbol in "
                             "replacement (must be 1-4 inclusive, found %d): "
                             "%s" % (len(symbol), symbol))
        for s,m in self.replacements:
            if s == symbol:
                raise ValueError("Must use a unique symbol for every isotopic "
                                 "replacement: %s" % s)
        mass = float(mass)
        if mass < 0.0 or mass > 500.0:
            raise ValueError("Unexpected isotopic replacement mass: %f" % mass)
        self.replacements.append((symbol,mass))

# read in atomic weight data
elements = []

root = os.path.split(os.path.abspath(__file__))[0]

with open(root + "/weights.dat", "r") as _weights_file:
    weights_lines = _weights_file.readlines()
for line in weights_lines:
    # ignore comments and blank lines
    line = line.strip()
    if len(line) == 0 or line[0] == "#":
        continue
    line = line.split("#",1)[0]

    # parse (the shipped weights.dat is well-formed, so the guards below are
    # defensive against a corrupted data file)
    fields = line.split(",")
    if len(fields) < 4:  # pragma: no cover
        raise ValueError("Not enough data on this line of weights.dat: %s"
                         % line)
    element = Element(*fields[0:4])
    if (len(fields)-4) % 2 != 0:  # pragma: no cover
        raise ValueError("Unexpected number of isotopic replacement fields on "
                         "this line of weights.dat (the number of fields after "
                         "the first four must be a multiple of 2, found %d): %s"
                         % (len(fields)-4, line))
    if (len(fields) > 4):
        for i in range(4, len(fields), 2):
            element.add_replacement(fields[i], fields[i+1])
    elements.append(element)
logger.info("Read atomic weight data for %d elements.", len(elements))

# map from atomic number to default masses
DEFAULT_MASSES = { e.atomic_number : e.default_mass for e in elements }

# map from valid isotopic replacements to masses
REPLACEMENTS = {}
for e in elements:
    for replacement,mass in e.replacements:
        REPLACEMENTS[replacement] = mass

# map from isotopic replacements to atomic numbers
REPLACEMENTS_Z = {}
for e in elements:
    for replacement,mass in e.replacements:
        REPLACEMENTS_Z[replacement]=e.atomic_number


def replacement_mass(replacement):
    """Resolve an isotopic replacement to a mass (amu).

    A number is used directly as the mass (for unusual/non-standard isotopes);
    a string is looked up among the standard isotope labels in weights.dat
    (e.g. "13C", "2D", "17O"). Raises ValueError for anything else.
    """
    if isinstance(replacement, bool):
        raise ValueError("invalid isotopic replacement: %r" % replacement)
    if isinstance(replacement, (int, float)):
        if replacement <= 0.0:
            raise ValueError("isotopic replacement mass must be positive: %r" % replacement)
        return float(replacement)
    if replacement in REPLACEMENTS:
        return REPLACEMENTS[replacement]
    raise ValueError("invalid isotopic replacement: %r" % replacement)


# threshold to separate linear molecules from non-linear molecules
LINEARITY_THRESHOLD = 1e-06
DROP_NUM_LINEAR = 5
# DROP_NUM_NONLINEAR = 6

# canonical element symbol -> atomic number map (the full periodic table; the
# masses in weights.dat only cover the subset PyQuiver ships isotope data for)
SYMBOL_TO_Z = {
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
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
}
