# This file reads 
import sys
import re
from constants import REPLACEMENTS, REPLACEMENTS_Z
from quiver import System

# Reads PyQuiver .config files.
class Config(object):
    def __init__(self,filename):
        expected_fields = "scaling temperature reference_isotopologue".split(" ")
        config = { i : None for i in expected_fields }
        config["filename"] = filename

        # a list of isotopologues
        # each entry is a list of tuples
        # each tuple is (from_atom_number, to_atom_number, replacement_isotope)
        # this format allows for multiple replacements in one isotopologue
        isotopologues = []
        
        # read file
        for line in open(filename, "r"):
            # ignore comments and blank lines
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue
            line = line.split("#",1)[0]

            # read space-delimited data
            fields = line.encode("ascii","ignore").split(" ")
            
            # for backwards compatibility with quiver
            if fields[0] == "isotopomer":
                fields[0] = "isotopologue"
            elif fields[0] == "reference_isotopomer":
                fields[0] = "reference_isotopologue"

            # parse
            if fields[0] == "isotopologue":
                # parse isotopologues, checking for sanity
                if len(fields) != 5:
                    raise ValueError("unexpected number of fields for isotopologue in config file:\n%s" % line)
                isotopologue_number, from_atom_number, to_atom_number, replacement = int(fields[1]), int(fields[2]), int(fields[3]), str(fields[4])
                if from_atom_number < 1 or to_atom_number < 1:
                    raise ValueError("check atom numbers:\n%s" % line)
                if replacement not in REPLACEMENTS.keys():
                    raise ValueError("invalid isotopic replacement:\n%s" % line)

                # allows for the fact that isotopologues can make multiple replacements
                if isotopologue_number == len(isotopologues) + 1:
                    # this is a new isotopologue
                    isotopologues.append([(from_atom_number, to_atom_number, replacement)])
                elif isotopologue_number == len(isotopologues):
                    # this is an additional replacement
                    isotopologues[-1].append((from_atom_number, to_atom_number, replacement))
                else:
                    # this is not allowed
                    raise ValueError("isotopologues must be numbered in ascending order:\n%s" % line)
            elif len(fields) == 2:
                # read regular configuration fields that have only one entry
                fields = [ str(i) for i in fields ]
                if fields[0] not in expected_fields:
                    raise ValueError("unexpected config field:\n%s" % line)
                config[fields[0]] = fields[1]
            else:
                raise ValueError("unexpected number of fields in config file:\n%s" % line)

        # ensure we have all the fields we are supposed to
        for k,v in config.iteritems():
            if v is None:
                raise ValueError("missing config file field: %s" % k)
        if len(isotopologues) == 0:
            raise ValueError("must specify at least one isotopologue")

        # check some of the other fields
        config["temperature"] = float(config["temperature"])
        if config["temperature"] < 0.0:
            raise ValueError("check temperature")
        config["scaling"] = float(config["scaling"])
        if config["scaling"] < 0.5 or config["scaling"] > 1.5:
            raise ValueError("check scaling factor")
        config["reference_isotopologue"] = int(config["reference_isotopologue"])
        if config["reference_isotopologue"] < 0 or config["reference_isotopologue"] > len(isotopologues):
            raise ValueError("check reference isotopologue is valid")
 
        # store all the information in the object dictionary
        config["isotopologues"] = isotopologues
        self.__dict__ = config

    # checks if this config file is compatible with a pair of ground and transition state systems
    def check(self, gs, ts, verbose=False):
        if not isinstance(gs, System):
            raise ValueError("gs is not a System")
        if not isinstance(ts, System):
            raise ValueError("ts is not a System")
        n_gs_atoms = len(gs.atomic_numbers)
        n_ts_atoms = len(ts.atomic_numbers)
        
        # check that the isotopic replacements make sense
        isotopologues = self.isotopologues
        for i,isotopologue in enumerate(isotopologues):
            isotopologue_number = i+1
            for j in range(len(isotopologue)):
                # this replacement changes gs atom number from_atom
                # and ts atom number to_atom
                # to replacement (like "2D")
                from_atom = isotopologue[j][0]
                to_atom = isotopologue[j][1]
                replacement = isotopologue[j][2]
                replacement_string = "%d %d %d %s" % (i+1, from_atom, to_atom, replacement)

                # get the atomic numbers of from_atom and to_atom
                # must subtract one to convert from atom numbers to indices
                from_atomZ = gs.atomic_numbers[from_atom-1]
                to_atomZ = gs.atomic_numbers[to_atom-1]
                replacementZ = REPLACEMENTS_Z[replacement]
                if from_atomZ != replacementZ:
                    raise ValueError("gs atomic number and replacement atomic number do not match for\n" % replacement_string)
                if to_atomZ != replacementZ:
                    raise ValueError("ts atomic number and replacement atomic number do not match for\n" % replacement_string)
                if (from_atomZ != to_atomZ):
                    raise ValueError("gs and ts atomic number do not match for\n" % replacement_string)

        if verbose:
            print "Config file %s makes sense with gs file %s and ts file %s.\n" % (self.filename, gs.filename, ts.filename)

    # convert to human-readable format
    def __str__(self):
        to_string = "Config file: %s\nTemperature: %.1f\nScaling: %.3f\nReference Isotopologue: %d\n" % (self.filename, self.temperature, self.scaling, self.reference_isotopologue)
        for i,isotopologue in enumerate(self.isotopologues):
            for j in range(len(isotopologue)):
                to_string += "   Isotopologue %2d, replacement %1d: replace gs atom %3d and ts atom %3d with %3s\n" % (i+1, j+1, isotopologue[j][0], isotopologue[j][1], isotopologue[j][2])
        return to_string[:-1]

#print config.__dict__

gs = System("../test/claisen_gs.out")
ts = System("../test/claisen_ts.out")
config = Config("test.config")
print config
config.check(gs, ts, verbose=True)

