import re

import numpy  as np
import pandas as pd

root = ""

def _g09_triangle_serial(row,col):
    if col > row:
        return _g09_triangle_serial(col,row)
    triangle = lambda n: n*(n+1)/2
    return triangle(row) + col

class Isotopologue(object):
    def __init__(self, system, masses):
        self.iitensor = self.calculate_inertia_tensor(masses)
        self.mw_hessian = self.calculate_mw_hessian(masses)
    
    def calculate_inertia_tensor(self, masses):
        # calculate cartesian center of mass to find intertia tensor relative to com
        rcm = np.zeros(3)
        total_mass = 0
        for i in xrange(0, number_of_atoms):
            for e in xrange(0,3):
                total_mass += masses[i]
                rcm[e] += positions[i][e] * masses[i]
        for e in xrange(0,3):
            rcm[e] = rcm[e] / total_mass

        # calculate cartesian moment of inertia tensor
        iitensor = np.zeros(shape=(3,3))
        for e1 in xrange(0,3):
            for e2 in xrange(0,3):
                for i in xrange(0, number_of_atoms):
                    if e1 == e2:
                        iitensor[e1,e2] += masses[i] * ((positions[i,(e1+1)%3]-rcm[(e1+1)%3])**2
                                                      + (positions[i,(e1+2)%3]-rcm[(e1+2)%3])**2)
                    else:
                        iitensor[e1,e2] += -1 * masses[i] * (positions[i,e1]-rcm[e1]) * (positions[i,e2]-rcm[e2])

        return iitensor

    def calculate_mw_hessian(self):
        pass
    
    def calculate_rpfr(self):
        pass


class System(object):
    def __init__(self, outfile, style="g09"):
        print outfile
        with open(outfile, 'r') as f:
            out_data = f.read()
            if style == "g09":
                # read in the number of atoms
                m = re.search("NAtoms\= + ([0-9]+)", out_data)
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
                    masses[center_number] = raw_geom_line[1]
                    for e in xrange(0,3):
                        positions[center_number][e] = raw_geom_line[3+e]
                                                

                # units = hartrees/bohr^2 ?
                hessian = self._parse_g09_hessian(out_data)

        #copy fields
        self.hessian = hessian

    def _parse_g09_hessian(self, data):
        #    number_of_atoms = 
        m = re.search("NImag\=(.+?)\@", data, re.DOTALL)
        if m:
            raw_archive = re.sub('[\s+]', '', m.group(0))
        else:
            raise AttributeError("No frequency job detected.")

        raw_fcm = raw_archive.split('\\')[2].split(',')
        fcm = np.zeros(shape=(3*self.number_of_atoms, 3*self.number_of_atoms))
        for i in xrange(0, 3*self.number_of_atoms):
            for j in xrange(0, 3*self.number_of_atoms):
                fcm[i,j] = raw_fcm[_g09_triangle_serial(i,j)]
        return fcm

gs = System("../test/claisen_gs.out")

