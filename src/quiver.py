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
                m = re.search("NAtoms\= + ([0-9]+)", out_data)
                if m:
                    number_of_atoms = int(m.group(1))
                else:
                    raise AttributeError("Number of atoms not detected.")
                self.number_of_atoms = number_of_atoms
                
                #m = re.match("Forces \(Hartrees/Bohr\)(.+)Cartesian Forces")
                for match in re.finditer("Forces \(Hartrees/Bohr\)(.+)Cartesian Forces", out_data, re.DOTALL):
                    print match.group(0)

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

