import numpy  as np
import pandas as pd
import re

from constants import DEFAULT_MASSES, PHYSICAL_CONSTANTS
from config import Config

def proj(u,v):
    # project u onto v
    return np.inner(v,u)/np.inner(u,u) * u

def normalize(v):
    norm=np.linalg.norm(v)
    if norm < 1.0E-5: 
       raise ValueError
    return v/norm

def test_orthogonality(vectors):
    mat = np.zeros(shape=(len(vectors),len(vectors)))
    for i,u in enumerate(vectors):
        for j,v in enumerate(vectors):
            inner = np.inner(u,v)
            if inner != 0.0 or inner != 1.0:
                mat[i][j] = inner
    print "Orthogonality:"
    for l in mat:
        print l

def schmidt(seed_vectors, rest_vectors, dimension):
    vectors = list(seed_vectors)
    test_vectors = list(rest_vectors)
    while len(vectors) < dimension:
        try:
            test_vector = test_vectors.pop()
            for v in vectors:
                test_vector -= proj(v, test_vector)
            try:
                vectors.append(normalize(test_vector))
            except ValueError:
                pass

        except IndexError:
            raise ValueError("Could not fill the appropriate dimensional space")            
    if len(vectors) < dimension:
        raise ValueError("Could not fill the appropriate dimensional space")
    else:
        return vectors

def _g09_triangle_serial(row,col):
    if col > row:
        return _g09_triangle_serial(col,row)
    triangle = lambda n: n*(n+1)/2
    return triangle(row) + col

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
        

        #self.frequencies = self.calculate_frequencies()
    
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
        

    def calculate_internal_hessian(self, masses):
        print "Masses:", masses
        print "RCM:", self.rcm
        print "Positions:", self.system.positions
        print "RCM Positions:", self.rcm_positions
        print "IITensor:", self.iitensor

        vectors = []
        for e in xrange(0,3):
            v = np.zeros(3*self.number_of_atoms)
            for i in xrange(0, self.number_of_atoms):
                v[3*i:3*i+3] = np.array([1.0 if x == e else 0.0 for x in xrange(0,3)]) * np.sqrt(masses[i])
            vectors.append(v)
        

        # order concerns?
        w,v = np.linalg.eig(self.iitensor)
        # which one?
        print "Eigenvalues:", w
        print "\"Eigenvectors\":", v
        x = np.column_stack(v)
        print "Stacked X matrix:", x

        mat = np.zeros(shape=(len(x),len(x)))
        
        bad_indices = []
        bad_vectors = []
        for i1 in xrange(len(w)):
            for i2 in xrange(i1+1, len(w)):
                if np.abs(w[i1]-w[i2]) < 1.0E-10:
                    bad_indices.extend([i1,i2])
                    bad_vectors.extend([v[:,i1], v[:,i2]])
        
        if bad_indices:
            print "Degenercies in eigenvalues. Assuming principal axes are not orthonormal and applying schmidt process."

        principal_axes = [None for i in w]
        for i in xrange(len(w)):
            if i not in bad_indices:
                principal_axes[i] = (v[:,i])
        

        adjusted_vectors = schmidt([], bad_vectors, len(bad_vectors))
        for i,v_ in enumerate(adjusted_vectors):
            principal_axes[bad_indices[i]] = v_
        
        print vectors
        print "Testing orthogonality of principal axes:"
        test_orthogonality(principal_axes)

        print principal_axes
        x = np.matrix.transpose(np.column_stack(principal_axes))
        print x

        # from the gaussian document: d=0 => making d4, d=1 => making d5 etc. j=j
        for d in xrange(0,3):
            v = np.zeros(3*self.number_of_atoms)
            axis = [1.0 if k == d else 0.0 for k in xrange(0,3)]
            for i in xrange(0, self.number_of_atoms):
                v[3*i:3*i+3] = np.cross(np.dot(x,self.rcm_positions[i]), axis)
            vectors.append(v)
            
        normalized_vectors = []
        zero_vectors = []
        for v in vectors:
            try:
                normalized_vectors.append(normalize(v))
            except ValueError:
                zero_vectors.append(v)
        print "Number of Zero Vectors:", len(zero_vectors)
        print "D Vectors:"
        for nv in normalized_vectors:
            print nv
        
        test_orthogonality(normalized_vectors)

        standard_basis = [np.array([1.0 if x == i else 0.0 for x in xrange(0,3*self.number_of_atoms)]) for i in xrange(0,3*self.number_of_atoms)]
        print normalized_vectors
        normalized_vectors = schmidt(normalized_vectors, standard_basis, 3*self.number_of_atoms)

        # costly step
        #print len(zero_vectors)
        d_matrix = np.matrix(np.column_stack(zero_vectors + normalized_vectors))
        
        # conversion factor to take hartree/(bohr^2 * amu) to units 1/s^2
        conv_factor = PHYSICAL_CONSTANTS['Eh']/(PHYSICAL_CONSTANTS['a0']**2 * PHYSICAL_CONSTANTS['amu'])
        
        # calculate the internal hessian in appropriate units
        int_hessian = np.dot(np.matrix.transpose(d_matrix), np.dot(self.mw_hessian, d_matrix)) * conv_factor
        return int_hessian

    def calculate_frequencies(self, threshold, method="mass weighted hessian"):
        # short circuit if frequencies have already been calculated
        print self.frequencies
        if self.frequencies is not None:
            return self.frequencies

        if method == "projected hessian":
            # need to detect if linear!!! TODO
            projected_hessian = self.int_hessian[np.ix_(range(5,3*self.number_of_atoms),range(5,3*self.number_of_atoms))]
            #projected_hessian = int_hessian[np.ix_(range(0,len(zero_vectors)) + range(6,3*self.number_of_atoms),range(0,len(zero_vectors)) + range(6,3*self.number_of_atoms))]

            #np.savetxt("int.csv", int_hessian, delimiter=",")
            #np.savetxt("proj.csv", projected_hessian, delimiter=",")

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
            return frequencies
        if method == "mass weighted hessian":
            imaginary_freqs = []
            small_freqs = []
            freqs = []

            v,w = np.linalg.eig(self.mw_hessian)
            frequencies = []
            for lam in v:
                print lam
                freq = np.sqrt(np.abs(lam))/(2*np.pi*PHYSICAL_CONSTANTS['c'])
                if np.linalg.norm(lam) < threshold:
                    small_freqs.append(freq)
                elif lam < 0: 
                    imaginary_freqs.append(-1 * freq)
                else:
                    freqs.append(freq)
                        
            frequencies = np.array(frequencies)
            frequencies.sort()
            self.frequencies = (small_freqs, imaginary_freqs, freqs)
            return self.frequencies


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

if __name__ == "__main__":
    #gs = System("../test/co2.out")
    #gsiso = Isotopologue(gs, gs.masses)
    from kie import KIE_Calculation
    KIE_Calculation("test.config", "../test/claisen_gs.out", "../test/claisen_ts.out")
#gsiso = Isotopologue(gs, gs.masses)
