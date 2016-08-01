import re

import numpy  as np
import pandas as pd

root = ""

def _g09_triangle_serial(row,col):
    if col > row:
        return _g09_triangle_serial(col,row)
    triangle = lambda n: n*(n+1)/2
    return triangle(row) + col

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
                atomic_number = [0 for i in xrange(number_of_atoms)]
                position = [None for i in xrange(number_of_atoms)]
                mass = np.array([0.0 for i in xrange(number_of_atoms)])
                for m in re.finditer("Standard orientation:(.+?)Rotational constants \(GHZ\)", out_data, re.DOTALL):
                    pass
                
                for l in m.group(1).split('\n')[5:-2]:
                    raw_geom_line = filter(None, l.split(' '))
                    center_number = int(raw_geom_line[0]) - 1
                    atomic_number[center_number] = int(raw_geom_line[1])
                    mass[center_number] = raw_geom_line[1]
                    position[center_number] = np.array(raw_geom_line[3:])

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
    
    def calculate_rpfr(self):
        pass


gs = System("../test/claisen_gs.out")

