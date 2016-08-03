import re

import numpy  as np
import pandas as pd
from constants import DEFAULT_MASSES
from constants import PHYSICAL_CONSTANTS

root = ""

def _g09_triangle_serial(row,col):
    if col > row:
        return _g09_triangle_serial(col,row)
    triangle = lambda n: n*(n+1)/2
    return triangle(row) + col

class IsotopolgueFactory(object):
    pass

class Isotopologue(object):
    def __init__(self, system, masses):
        self.system = system
        self.masses = masses

        # create vector of masses with each entry repeated three times for convenience
        masses3_list = []
        for m in masses:
            masses3_list.extend([m, m, m])
        self.masses3 = np.array(masses3_list)

        self.number_of_atoms = system.number_of_atoms
        self.rcm, self.rcm_positions, self.iitensor = self.calculate_inertia_tensor(masses, system.positions)
        #self.rcm, self.iitensor = self.calculate_inertia_tensor(masses, system.positions_angstrom)
        self.mw_hessian = self.calculate_mw_hessian(self.masses3)

        self.calculate_internal_hessian(masses)

    
    def calculate_inertia_tensor(self, masses, positions):
        # calculate cartesian center of mass to find intertia tensor relative to center of mass
        rcm = np.zeros(3)
        total_mass = 0
        #import pdb; pdb.set_trace()
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

    def calculate_frequencies(self):
        pass

    def calculate_mw_hessian(self, masses3):
        hessian = self.system.hessian
        mw_hessian = np.zeros_like(hessian)

        for i in xrange(0, 3*self.number_of_atoms):
            for j in xrange(0, 3*self.number_of_atoms):
                mw_hessian[i,j] = hessian[i,j] / np.sqrt( masses3[i] * masses3[j] )

        return mw_hessian
                

    def calculate_internal_hessian(self, masses):
        vectors = []
        for e in xrange(0,3):
            v = np.zeros(3*self.number_of_atoms)
            for i in xrange(0, self.number_of_atoms):
                v[3*i:3*i+3] = np.array([1.0 if x == e else 0.0 for x in xrange(0,3)]) * np.sqrt(masses[i])
            vectors.append(v)
        

        # order concerns?
        _,w = np.linalg.eig(self.iitensor)

        # which one?
        #x = np.vstack(w)
        x = np.column_stack(w)
        '''
        for k in xrange(0,3):
            v = []
            for i in xrange(0, self.number_of_atoms):
                di = np
                v.extend(di)
            vectors.append(np.array(v))
        '''

        # from the gaussian document: d=0 => making d4, d=1 => making d5 etc. j=j
        for d in xrange(0,3):
            v = np.zeros(3*self.number_of_atoms)
            for i in xrange(0, self.number_of_atoms):
                p = np.zeros(3)
                for j in xrange(0,3):
                    p[j] = np.inner(x[:, j], self.rcm_positions[i])

                for j in xrange(0,3):
                    v[3*i+j] = (p[(d+1)%3] * x[j][(d+2)%3] - p[(d+2)%3] * x[j][(d+1)%3])/np.sqrt(masses[i])
            vectors.append(v)

        def normalize(v):
            norm=np.linalg.norm(v)
            if norm==0: 
               raise ValueError
            return v/norm

        normalized_vectors = []
        zero_vectors = []
        for v in vectors:
            try:
                normalized_vectors.append(normalize(v))
            except ValueError:
                zero_vectors.append(v)

        for u in normalized_vectors:
            for v in normalized_vectors:
                pass
                #print np.inner(u,v)

        #print normalized_vectors

        def proj(u,v):
            # project u onto v
            #print np.inner(v,u)/np.inner(u,u)
            return np.inner(v,u)/np.inner(u,u) * u

        standard_basis = [np.array([1.0 if x == i else 0.0 for x in xrange(0,3*self.number_of_atoms)]) for i in xrange(0,3*self.number_of_atoms)]
        while len(normalized_vectors) + len(zero_vectors) < 3*self.number_of_atoms:
            test_vector = standard_basis.pop()
            for v in normalized_vectors:
                test_vector -= proj(v, test_vector)
            try:
                normalized_vectors.append(normalize(test_vector))
            except ValueError:
                pass
        
        # tiny residuals left over but otherwise good
        '''
        for u in normalized_vectors:
            for v in normalized_vectors:
                print np.inner(u,v)
        '''

        # costly step
        print len(zero_vectors)
        d_matrix = np.matrix(np.column_stack(zero_vectors + normalized_vectors))
        # conversion factor to take hartree/(bohr^2 * amu) to units 1/s^2
        conv_factor = PHYSICAL_CONSTANTS['Eh']/(PHYSICAL_CONSTANTS['a0']**2 * PHYSICAL_CONSTANTS['amu'])
        
        # calculate the internal hessian in appropriate units
        int_hessian = np.dot(np.matrix.transpose(d_matrix), np.dot(self.mw_hessian, d_matrix)) * conv_factor

        # need to detect if linear!!! TODO
        projected_hessian = int_hessian[np.ix_(range(6,3*self.number_of_atoms),range(6,3*self.number_of_atoms))]
        #projected_hessian = int_hessian[np.ix_(range(0,len(zero_vectors)) + range(6,3*self.number_of_atoms),range(0,len(zero_vectors)) + range(6,3*self.number_of_atoms))]

        np.savetxt("int.csv", int_hessian, delimiter=",")
        np.savetxt("proj.csv", projected_hessian, delimiter=",")

        v,w = np.linalg.eig(projected_hessian)
        #v,w = np.linalg.eig(self.mw_hessian*conv_factor)
        #v,w = np.linalg.eig(int_hessian)
        # retrieve frequencies in units 1/cm
        frequencies = []
        for lam in v:
            imag_flag = 1
            if lam < 0:
                imag_flag = -1
            lam = np.abs(lam)
            frequencies.append(imag_flag * np.sqrt(lam)/(2*np.pi*PHYSICAL_CONSTANTS['c']))
        frequencies = np.array(frequencies)
        frequencies.sort()
        print frequencies
        return int_hessian
        
    def calculate_rpfr(self):
        pass


class System(object):
    def __init__(self, outfile, style="g09"):
        print "Reading data from %s..." % outfile
        self.filename = outfile
        with open(outfile, 'r') as f:
            out_data = f.read()
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
                    masses[center_number] = DEFAULT_MASSES[atomic_numbers[center_number]]
                    for e in xrange(0,3):
                        positions[center_number][e] = raw_geom_line[3+e]
                
                # units = hartrees/bohr^2 ?
                hessian = self._parse_g09_hessian(out_data)

        #copy fields
        self.hessian = hessian

        self.positions_angstrom = positions
        self.positions = positions/PHYSICAL_CONSTANTS['atb']
        #print self.positions_angstrom[0]
        #print self.positions[0]

        self.masses = masses
        self.atomic_numbers = atomic_numbers

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
                fcm[i,j] = raw_fcm[_g09_triangle_serial(i,j)]
        return fcm

gs = System("../test/claisen_gs.out")
gsiso = Isotopologue(gs, gs.masses)
