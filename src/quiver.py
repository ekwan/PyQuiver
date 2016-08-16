import argparse
import sys
import re
import os

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
        self.mw_hessian = self.calculate_mw_hessian(self.masses3)

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
            #v,w = np.linalg.eigh(self.mw_hessian*conv_factor)
            v = np.linalg.eigvalsh(self.mw_hessian*conv_factor)
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
                lines = out_data.split("\n")
                try:
                    number_of_atoms = int(lines[0])
                except ValueError:
                    raise ValueError("first line must contain integer number of atoms.")
                self.number_of_atoms = number_of_atoms

                atomic_numbers = [0 for i in xrange(number_of_atoms)]
                positions = np.zeros(shape=(number_of_atoms,3))
                
                for l in lines[1:number_of_atoms+1]:
                    fields = l.split(',')
                    try:
                        center_number, atomic_number, x, y, z = fields
                    except ValueError:
                        raise ValueError("the following line in the geometry did not have the appropriate number of fields: {0}".format(l))
                    center_number = int(center_number)
                    atomic_number = int(atomic_number)

                    atomic_numbers[center_number] = atomic_number
                    positions[center_number][0] = x
                    positions[center_number][1] = y
                    positions[center_number][2] = z
                
                fcm_fields = lines[number_of_atoms+1].split(',')
                hessian = self._parse_serial_lower_hessian(fcm_fields)
                    
            if style == "g09":
                # read in the number of atoms
                m = re.search("NAtoms\= +([0-9]+)", out_data)
                if m:
                    number_of_atoms = int(m.group(1))
                else:
                    raise AttributeError("Number of atoms not detected.")
                m = None
                self.number_of_atoms = number_of_atoms
                # read in the last geometry (assumed cartesian coordinates)
                atomic_numbers = [0 for i in xrange(number_of_atoms)]
                positions = np.zeros(shape=(number_of_atoms,3))
                for m in re.finditer("Standard orientation(.+?)Rotational constants \(GHZ\)", out_data, re.DOTALL):
                    pass

                if not m:
                    raise AttributeError("Geometry table not detected.")

                def valid_geom_line_p(split_line):
                    if len(split_line) == 6:
                        try:
                            int(split_line[0])
                            int(split_line[1])
                            int(split_line[2])
                            return True
                        except ValueError:
                            return False
                    return False

                for l in m.group(1).split('\n'):
                    raw_geom_line = l.split()
                    raw_geom_line = filter(None, l.split(' '))
                    #print raw_geom_line
                    if valid_geom_line_p(raw_geom_line):
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

        self.atomic_numbers = atomic_numbers

    def _lower_triangle_serial_triangle(self,row,col):
        if col > row:
            return self._lower_triangle_serial_triangle(col,row)
        triangle = lambda n: n*(n+1)/2
        return triangle(row) + col

    def _parse_g09_hessian(self, data):
        raw_archive = re.findall(".*l9999.exe(.+?)\@", data, re.DOTALL)[-1]
        raw_archive = re.sub('[\s+]', '', raw_archive) + "@"
        m = re.search("NImag\=(.+?)\@", raw_archive, re.DOTALL)
        if m:
            pass
        else:
            raise AttributeError("No frequency job detected.")

        raw_fcm = m.group(0).split('\\')[2].split(',')
        fcm = self._parse_serial_lower_hessian(raw_fcm)
        return fcm

    def _parse_serial_lower_hessian(self, fields):
        fcm = np.zeros(shape=(3*self.number_of_atoms, 3*self.number_of_atoms))
        for i in xrange(0, 3*self.number_of_atoms):
            for j in xrange(0, 3*self.number_of_atoms):
                fcm[i,j] = fields[self._lower_triangle_serial_triangle(i,j)]
        return fcm

    def _make_serial_hessian(self):
        serial = ""
        print self.hessian
        print self.hessian[1]
        for i in xrange(self.number_of_atoms*3):
            for j in xrange(0,i+1):
                #print self.hessian[i,j]
                serial += str(self.hessian[i,j]) + ","
        return serial

    def _make_serial_geometry(self):
        serial = ""
        for i in xrange(self.number_of_atoms):
            serial += "{0},{1},{2},{3},{4}\n".format(i, self.atomic_numbers[i], self.positions_angstrom[i,0], self.positions_angstrom[i,1], self.positions_angstrom[i,2])
        return serial


    def dump_pyquiver_input_file(self, extension=".qin"):
        path = os.path.splitext(self.filename)[0] + extension
        serial = str(self.number_of_atoms) + "\n"
        serial += self._make_serial_geometry()
        serial += self._make_serial_hessian()

        with open(path, 'w') as f:
            f.write(serial)
        
        return serial

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A program that calculates KIEs and EIEs based on a ground and transition state file.")
    parser.add_argument('-v', '--verbose', dest="debug", help='when the verbose flag is set debug information is printed', action='store_true')
    parser.add_argument('-s', '--style', dest="style", default='g09', help='style of input files')
    parser.add_argument('config', help='configuration file path')
    parser.add_argument('gs', help='ground state file path')
    parser.add_argument('ts', help='transition state file path')

    args = parser.parse_args()

    from kie import KIE_Calculation
    calc = KIE_Calculation(args.config, args.gs, args.ts, style=args.style, debug=args.debug)
    print calc
