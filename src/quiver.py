import sys
import re

import numpy  as np

from utility import proj, normalize, test_orthogonality, schmidt
from constants import DEFAULT_MASSES, PHYSICAL_CONSTANTS
from config import Config

# represents a geometric arrangement of atoms with specific masses
class Isotopologue(object):
    def __init__(self, id_, system, masses):
        self.name = id_
        self.system = system
        self.masses = masses
        self.frequencies = None
        # create vector of masses with each entry repeated three times for convenience
        masses3_list = []
        for m in masses:
            masses3_list.extend([m, m, m])
        self.masses3 = np.array(masses3_list)

        self.number_of_atoms = system.number_of_atoms
        self.rcm, self.rcm_positions, self.iitensor = self.calculate_inertia_tensor(masses, system.positions)
        self.mw_hessian = self.calculate_mw_hessian(self.masses3)
        return
        self.int_hessian = self.calculate_internal_hessian(masses)
        
    def calculate_inertia_tensor(self, masses, positions):
        # calculate cartesian center of mass to find intertia tensor relative to center of mass
        rcm = np.zeros(3)
        total_mass = 0
        for i in xrange(0, self.number_of_atoms):
            total_mass += masses[i]
            rcm += positions[i] * masses[i]

        rcm = rcm / total_mass

        # calculate cartesian moment of inertia tensor
        iitensor = np.zeros(shape=(3,3))
        # center to rcm and convert to atomic units
        rcm_positions = positions - rcm

        for e1 in xrange(0,3):
            for e2 in xrange(0,3):
                for i in xrange(0, self.number_of_atoms):
                    if e1 == e2:
                        iitensor[e1,e2] += masses[i] * ((rcm_positions[i,(e1+1)%3])**2
                                                      + (rcm_positions[i,(e1+2)%3])**2)
                    else:
                        iitensor[e1,e2] += -1 * masses[i] * (rcm_positions[i,e1]) * (rcm_positions[i,e2])

        return rcm, rcm_positions, iitensor

    def calculate_mw_hessian(self, masses3):
        hessian = self.system.hessian
        mw_hessian = np.zeros_like(hessian)

        for i in xrange(0, 3*self.number_of_atoms):
            for j in xrange(0, 3*self.number_of_atoms):
                mw_hessian[i,j] = hessian[i,j] / np.sqrt( masses3[i] * masses3[j] )

        return mw_hessian
        
    def calculate_frequencies(self, threshold, scaling=1.0, method="mass weighted hessian"):
        # short circuit if frequencies have already been calculated
        if self.frequencies is not None:
            return self.frequencies

        if method == "mass weighted hessian":
            imaginary_freqs = []
            small_freqs = []
            conv_factor = PHYSICAL_CONSTANTS['Eh']/(PHYSICAL_CONSTANTS['a0']**2 * PHYSICAL_CONSTANTS['amu'])
            v,w = np.linalg.eigh(self.mw_hessian*conv_factor)
            freqs = []

            for lam in v:
                freq = np.sqrt(np.abs(lam))/(2*np.pi*PHYSICAL_CONSTANTS['c'])*(scaling)
                if np.linalg.norm(freq) < threshold:
                    small_freqs.append(freq)
                elif lam < 0: 
                    imaginary_freqs.append(-1 * freq)
                else:
                    freqs.append(freq)

            imaginary_freqs.sort()
            freqs.sort()
            small_freqs.sort()
            if len(imaginary_freqs) > 1:
                print "WARNING: more than one imaginary frequency detected. Taking mode with largest norm."
                imaginary_freqs = [imaginary_freqs[-1]]
            self.frequencies = (small_freqs, imaginary_freqs, freqs)
            return self.frequencies


class System(object):
    def __init__(self, outfile, style="g09"):
        valid_styles = ["g09", "pyquiver"]
        if style not in valid_styles:
            raise ValueError("specified style, {0}, not supported".format(style))

        print "Reading data from {0}... with style {1}".format(outfile, style)
        self.filename = outfile
        with open(outfile, 'r') as f:
            out_data = f.read()
            if style == "pyquiver":
                pass

            if style == "g09":
                # read in the number of atoms
                m = re.search("NAtoms\= +([0-9]+)", out_data)
                if m:
                    number_of_atoms = int(m.group(1))
                else:
                    raise AttributeError("Number of atoms not detected.")
                self.number_of_atoms = number_of_atoms
                
                # read in the last geometry (assumed cartesian coordinates)
                atomic_numbers = [0 for i in xrange(number_of_atoms)]
                positions = np.zeros(shape=(number_of_atoms,3))
                masses = np.zeros(number_of_atoms)
                for m in re.finditer("Standard orientation:(.+?)Rotational constants \(GHZ\)", out_data, re.DOTALL):
                    pass
                
                for l in m.group(1).split('\n')[5:-2]:
                    raw_geom_line = filter(None, l.split(' '))
                    center_number = int(raw_geom_line[0]) - 1
                    atomic_numbers[center_number] = int(raw_geom_line[1])
                    for e in xrange(0,3):
                        positions[center_number][e] = raw_geom_line[3+e]
                
                # units = hartrees/bohr^2 ?
                hessian = self._parse_g09_hessian(out_data)

        #copy fields
        self.hessian = hessian

        self.positions_angstrom = positions
        self.positions = positions/PHYSICAL_CONSTANTS['atb']

        self.masses = masses
        self.atomic_numbers = atomic_numbers

    def _g09_triangle_serial(self,row,col):
        if col > row:
            return self._g09_triangle_serial(col,row)
        triangle = lambda n: n*(n+1)/2
        return triangle(row) + col

    def _parse_g09_hessian(self, data):
        m = re.search("NImag\=(.+?)\@", data, re.DOTALL)
        if m:
            raw_archive = re.sub('[\s+]', '', m.group(0))
        else:
            raise AttributeError("No frequency job detected.")

        raw_fcm = raw_archive.split('\\')[2].split(',')
        fcm = np.zeros(shape=(3*self.number_of_atoms, 3*self.number_of_atoms))
        for i in xrange(0, 3*self.number_of_atoms):
            for j in xrange(0, 3*self.number_of_atoms):
                fcm[i,j] = raw_fcm[self._g09_triangle_serial(i,j)]
        return fcm

if __name__ == "__main__":
    input_style = None
    print sys.argv
    if len(sys.argv) == 5:
        input_style = sys.argv.pop()
    if len(sys.argv) != 4:
        print "Usage: python quiver.py configuration_file ground_state_file transition_state_file"

    from kie import KIE_Calculation
    if input_style:
        KIE_Calculation(sys.argv[1], sys.argv[2], sys.argv[3], style=input_style)
    else:
        KIE_Calculation(sys.argv[1], sys.argv[2], sys.argv[3])
