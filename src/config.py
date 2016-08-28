# This file reads PyQuiver configuration files. 
import sys
import re
from constants import REPLACEMENTS, REPLACEMENTS_Z

# Reads PyQuiver .config files.
class Config(object):
    def __init__(self,filename):
        expected_fields = "scaling temperature mass_override_isotopologue reference_isotopologue frequency_threshold".split(" ")
        config = { i : None for i in expected_fields }
        config["filename"] = filename

        print "\nReading configuration from {0}".format(filename)

        # a list of isotopologues
        # each entry is a list of tuples
        # each tuple is (from_atom_number, to_atom_number, replacement_isotope)
        # this format allows for multiple replacements in one isotopologue
        isotopologues = {}
        
        # read file
        for line in open(filename, "r"):
            # ignore comments and blank lines
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue
            line = line.split("#", 1)[0]

            # read space-delimited data
            fields = filter(None, line.encode("ascii","ignore").split(" "))
            
            # for backwards compatibility with quiver
            if fields[0] == "isotopomer":
                fields[0] = "isotopologue"
            elif fields[0] == "reference_isotopomer":
                fields[0] = "reference_isotopologue"
            elif fields[0] == "mass_override_isotopomer":
                fields[0] = "mass_override_isotopologue"


            # parse
            if fields[0] == "isotopologue":
                # parse isotopologues, checking for sanity
                if len(fields) != 5:
                    raise ValueError("unexpected number of fields for isotopologue in config file:\n%s" % line)
                isotopologue_id, from_atom_number, to_atom_number, replacement = str(fields[1]), int(fields[2]), int(fields[3]), str(fields[4])
                if isotopologue_id == "default":
                    raise ValueError("name default is reserved.")
                if from_atom_number < 1 or to_atom_number < 1:
                    raise ValueError("check atom numbers:\n%s" % line)
                if replacement not in REPLACEMENTS.keys():
                    raise ValueError("invalid isotopic replacement:\n%s" % line)

                # allows for the fact that isotopologues can make multiple replacements
                try:
                    isotopologues[isotopologue_id].append((from_atom_number, to_atom_number, replacement))
                except KeyError:
                    isotopologues[isotopologue_id] = [(from_atom_number, to_atom_number, replacement)]

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

        try:
            config["reference_isotopologue"]
        except KeyError:
            raise ValueError("check reference isotopologue is valid")
        config["reference_isotopologue"] = str(config["reference_isotopologue"])

        try:
            config["mass_override_isotopologue"]
        except KeyError:
            raise ValueError("check reference isotopologue is valid")
        config["mass_override_isotopologue"] = str(config["mass_override_isotopologue"])

        config["frequency_threshold"] = float(config["frequency_threshold"])
        if config["frequency_threshold"] > 100.0:
            raise ValueError("frequency threshold is too high")
 
        # store all the information in the object dictionary
        config["isotopologues"] = isotopologues
        self.__dict__ = config

    # checks if this config file is compatible with a pair of ground and transition state systems
    def check(self, gs, ts, verbose=False):
        n_gs_atoms = len(gs.atomic_numbers)
        n_ts_atoms = len(ts.atomic_numbers)
        # check that the isotopic replacements make sense
        isotopologues = self.isotopologues
        for i,isotopologue in isotopologues.iteritems():
            for r in range(len(isotopologue)):
                # this replacement changes gs atom number from_atom
                # and ts atom number to_atom
                # to replacement (like "2D")
                from_atom = isotopologue[r][0]
                to_atom = isotopologue[r][1]
                replacement = isotopologue[r][2]
                replacement_string = "%s %d %d %s" % (i, from_atom, to_atom, replacement)
                # get the atomic numbers of from_atom and to_atom
                # must subtract one to convert from atom numbers to indices
                from_atomZ = gs.atomic_numbers[from_atom-1]
                to_atomZ = ts.atomic_numbers[to_atom-1]
                replacementZ = REPLACEMENTS_Z[replacement]
                if from_atomZ != replacementZ:
                    raise ValueError("gs atomic number and replacement atomic number do not match for {0}\n".format(replacement_string))
                if to_atomZ != replacementZ:
                    raise ValueError("ts atomic number and replacement atomic number do not match for {0}\n".format(replacement_string))
                if (from_atomZ != to_atomZ):
                    raise ValueError("gs and ts atomic number do not match for\n" % replacement_string)

        if verbose:
            print "Config file %s makes sense with gs file %s and ts file %s.\n" % (self.filename, gs.filename, ts.filename)

    # convert to human-readable format
    def __str__(self):
        to_string = "Config file: %s\nTemperature: %.1f K\nScaling: %.3f\nReference Isotopologue: %s\nFrequency threshold (cm-1): %d\n" % \
                    (self.filename, self.temperature, self.scaling, self.reference_isotopologue, self.frequency_threshold)

        keys = self.isotopologues.keys()
        if self.reference_isotopologue != "default" and self.reference_isotopologue != "none":
            try:
                keys.remove(self.reference_isotopologue)
            except:
                print "\nCould not find the following reference isotopologue: %s" % self.reference_isotopologue
                sys.exit(1)
        keys.sort()

        if self.reference_isotopologue == "default" or self.reference_isotopologue == "none":
            to_string += "   No reference isotopologue.\n"
        else:
            keys = [self.reference_isotopologue] + keys
        for i in keys:
            isotopologue = self.isotopologues[i]

            for j in range(len(isotopologue)):
                to_string += "   Isotopologue {0: >10s}, replacement {1: >2d}: replace gs atom {2: ^3d} and ts atom {3: ^3d} with {4: >3s}\n".format(i, j+1, isotopologue[j][0], isotopologue[j][1], isotopologue[j][2])

        return to_string[:-1]
