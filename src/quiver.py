import argparse
import sys
import re
import os
import pickle

import numpy as np

import settings
from utility import proj, normalize, test_orthogonality, schmidt
from constants import DEFAULT_MASSES, PHYSICAL_CONSTANTS, LINEARITY_THRESHOLD, DROP_NUM_LINEAR
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

    def __str__(self):
        returnString  = "Isotopologue: %s\n" % self.name
        returnString += " masses: %s" % self.masses.__str__()
        return returnString

    def dump_debug(self, name, obj):
        pass
        #f = open(name+'_'+(self.system.filename)+'.pickle', 'w')
        #ret = pickle.dump(obj, f)
        #f.close()

    def calculate_mw_hessian(self, masses3):
        hessian = self.system.hessian
        mw_hessian = np.zeros_like(hessian)

        mass_weights=[]

        for i in range(0, 3*self.number_of_atoms):
            for j in range(0, 3*self.number_of_atoms):
                mass_weights.append(1/(np.sqrt(masses3[i]*masses3[j])))
                mw_hessian[i,j] = hessian[i,j] / np.sqrt( masses3[i] * masses3[j] )
        if settings.DEBUG >= 3:
            self.dump_debug("mw", mass_weights)
            self.dump_debug("mw_hessian", mw_hessian)
        return mw_hessian

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
            freqs = [ np.copysign(np.sqrt(np.abs(freq)),freq) * constant for freq in v ]
            freqs.sort()

            imaginary_freqs = []
            small_freqs = []
            regular_freqs = []

            # detect imaginary frequencies
            for f in freqs:
                if f < -imag_threshold:
                    imaginary_freqs.append(f)

            if len(imaginary_freqs) > 1:
                print("WARNING: multiple imaginaries")

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
            if settings.DEBUG >= 3:
                self.dump_debug("freqs", self.frequencies)
            return self.frequencies
        else:
            raise ValueError("unknown frequency calculation type")


class System(object):
    def __init__(self, outfile, style="g09"):
        self.positions_angstrom = False
        self.positions = False
        self.atomic_numbers = False
        self.number_of_atoms = False
        self.hessian = False
        self.is_linear = True
        self.filename = outfile

        valid_styles = ["g09", "pyquiver", "orca"]
        if style not in valid_styles:
            raise ValueError("specified style, {0}, not supported".format(style))

        if settings.DEBUG >= 1:
            print("Reading data from {0}... with style {1}".format(outfile, style))

        # assumes snip files worked correctly
        if style == "g09" and not outfile.endswith(".snip") and not "Normal termination" in tail(outfile):
            raise ValueError("Gaussian job %s terminated in an error" % outfile)

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

                atomic_numbers = [0 for i in range(number_of_atoms)]
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
                # check that the verbose output option has been set
                verbose_flag_present = re.search(r" *#[pP] ", out_data)
                if verbose_flag_present == None and not outfile.lower().endswith(".snip"):
                    print()
                    print("Error: Gaussian output file %s" % outfile)
                    print("was not run with the verbose flag, so it does not contain enough information for")
                    print("PyQuiver to run.  Please re-run this calculation with a route card that starts with #p")
                    print()
                    sys.exit(1)

                # read in the number of atoms
                m = re.search("NAtoms\= +([0-9]+)", out_data)
                if m:
                    number_of_atoms = int(m.group(1))
                else:
                    raise AttributeError("Number of atoms not detected.")
                m = None
                self.number_of_atoms = number_of_atoms

                # read in the last geometry (assumed cartesian coordinates)
                atomic_numbers = [0 for i in range(number_of_atoms)]
                positions = np.zeros(shape=(number_of_atoms,3))

                # use standard orientation if possible
                for m in re.finditer("Standard orientation(.+?)Rotational constants \(GHZ\)", out_data, re.DOTALL):
                    pass

                # for input files with nosymm keyword, use input orientation
                if m is None:
                    for m in re.finditer("Input orientation(.+?)Distance matrix", out_data, re.DOTALL):
                        pass
                    if not m is None:
                        print("Couldn't find standard orientation so used input orientation instead.")

                if m is None:
                    for m in re.finditer("Input orientation(.+?)Rotational constants \(GHZ\)", out_data, re.DOTALL):
                        pass
                    if not m is None:
                        print("Couldn't find standard orientation so used input orientation instead.")

                # still couldn't find any geometries
                if m is None:
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
                    raw_geom_line = [_f for _f in l.split(' ') if _f]

                    if valid_geom_line_p(raw_geom_line):
                        center_number = int(raw_geom_line[0]) - 1
                        atomic_numbers[center_number] = int(raw_geom_line[1])
                        for e in range(0,3):
                            positions[center_number][e] = raw_geom_line[3+e]

                # units = hartrees/bohr^2
                hessian = self._parse_g09_hessian(out_data)

            elif style == "orca":
                from orca import parse_orca_output
                atomic_numbers, positions, hessian = parse_orca_output(out_data)
                self.number_of_atoms = len(atomic_numbers)


        # copy fields
        self.hessian = hessian

        self.positions_angstrom = positions
        # convert position in angstroms to needed bohr
        self.positions = positions/PHYSICAL_CONSTANTS['atb']

        # detect if molecule is linear
        # method: take every difference of two centers (that share center 0)
        #   calculate the dot product with other differences

        if self.number_of_atoms > 2:
            detected_indep = 0
            for i in range(1,len(positions)-1):
                for j in range(i+1, len(positions)):
                # compares how parallel the unit vectors are
                    diff0i = positions[0] - positions[i]
                    diff0j = positions[0] - positions[j]
                    u = diff0i/np.linalg.norm(diff0i)
                    v = diff0j/np.linalg.norm(diff0j)
                    # if the vectors are (probably) linearly indep, we break
                    if 1 - np.abs(np.inner(u,v)) > LINEARITY_THRESHOLD:
                        self.is_linear = False
                        detected_indep = 1
                        break

                if detected_indep:
                    break
        linear_string = "linear" if self.is_linear else "not linear"
        if settings.DEBUG >= 2:
            print("Molecule is %s." % linear_string)
        self.atomic_numbers = atomic_numbers

    def dump_debug(self, obj):
        f = open('freqs_'+self.name+'.json', 'w')
        out = json.dumps(self.frequencies)
        f.write(out)
        f.close()

    def _lower_triangle_serial_triangle(self,row,col):
        if col > row:
            return self._lower_triangle_serial_triangle(col,row)
        triangle = lambda n: n*(n+1)/2
        return int(triangle(row) + col)

    # search for the archive at the end of the file
    # then extract the force constant matrix
    def _parse_g09_hessian(self, data):
        #print("\n\nparsing\n\n")
        # regex for finding text between 1\1\GINC and \@
        # DOTALL means that . will match newlines
        # there are two capture groups here, which is why we have to use archive[0] later
        raw_archive = re.findall(r"1\\1\\GINC(.+?)\\(\s*)@", data, re.DOTALL)
        found_frequencies = False
        for archive in raw_archive:
            archive = re.sub('[\s+]', '', archive[0])
            #print(archive[:1000])
            #print("...")
            #print(archive[-1000:])
            #print("---")
            archive = re.search("NImag\=(.+?)$", archive, re.DOTALL)
            #print(archive)
            #print("*")
            #print()
            if archive:
                found_frequencies = True
                break
        if not found_frequencies:
            raise ValueError(f"No frequency job detected in {self.filename}.")

        raw_fcm = archive.group(0).split('\\')[2].split(',')
        #print(raw_fcm)
        self.raw_fcm = raw_fcm
        fcm = self._parse_serial_lower_hessian(raw_fcm)
        #print("\n\nsuccess\n\n")
        return fcm

    def _parse_serial_lower_hessian(self, fields):
        fcm = np.zeros(shape=(3*self.number_of_atoms, 3*self.number_of_atoms))
        for i in range(0, 3*self.number_of_atoms):
            for j in range(0, 3*self.number_of_atoms):
                fcm[i,j] = fields[self._lower_triangle_serial_triangle(i,j)]
        return fcm

    def _make_serial_hessian(self):
        serial = ""
        for i in range(self.number_of_atoms*3):
            for j in range(0,i+1):
                serial += str(self.hessian[i,j]) + ","
        return serial

    def _make_serial_geometry(self):
        serial = ""
        for i in range(self.number_of_atoms):
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
    parser = argparse.ArgumentParser(description="PyQuiver calculates KIEs and EIEs based on a ground and transition state file.")
    parser.add_argument('-v', '--verbose', dest="debug", help='when the verbose flag is set debug information is printed', action='count')
    parser.add_argument('-s', '--style', dest="style", default='g09', help='style of input files (g09, orca, or pyquiver)')
    parser.add_argument('config', help='configuration file path')
    parser.add_argument('gs', help='ground state file path')
    parser.add_argument('ts', help='transition state file path')

    args = parser.parse_args()
    if args.debug:
        settings.DEBUG = args.debug + 1

    from kie import KIE_Calculation
    calc = KIE_Calculation(args.config, args.gs, args.ts, style=args.style)
    print(calc)


def slugify(value):
    return "".join(x for x in value if x.isalnum())

def tail(filename):
    with open(filename) as f:
        content = f.readlines()
    return content[-1].strip()
