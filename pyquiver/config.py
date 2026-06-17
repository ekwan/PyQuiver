# This file reads PyQuiver configuration files.
import logging
from .constants import REPLACEMENTS_Z, replacement_mass
from collections import OrderedDict

logger = logging.getLogger("pyquiver")

EXPECTED_FIELDS = "scaling temperature mass_override_isotopologue reference_isotopologue imag_threshold frequency_threshold".split(" ")


# Reads PyQuiver .config files.
class Config(object):
    def __init__(self,filename):
        config = { i : None for i in EXPECTED_FIELDS }
        config["filename"] = filename

        logger.debug("Reading configuration from %s", filename)

        # a list of isotopologues
        # each entry is a list of tuples
        # each tuple is (from_atom_number, to_atom_number, replacement_isotope)
        # this format allows for multiple replacements in one isotopologue
        isotopologues = OrderedDict()
        expected_fields = EXPECTED_FIELDS

        # read file
        with open(filename, "r") as _config_file:
            config_lines = _config_file.readlines()
        for line in config_lines:
            # ignore comments and blank lines
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue
            line = line.split("#", 1)[0]

            # read space-delimited data
            fields = [_f for _f in line.split(" ") if _f]

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
                isotopologue_id, from_atom_number, to_atom_number = str(fields[1]), int(fields[2]), int(fields[3])
                # a replacement is a standard isotope label (e.g. 13C) or a
                # bare number for an unusual mass (e.g. 5000.0)
                try:
                    replacement = float(fields[4])
                except ValueError:
                    replacement = str(fields[4])
                if isotopologue_id == "default":
                    raise ValueError("name default is reserved.")
                if from_atom_number < 1 or to_atom_number < 1:
                    raise ValueError("check atom numbers:\n%s" % line)
                try:
                    replacement_mass(replacement)   # validates label or numeric mass
                except ValueError:
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

        self._finalize(config, isotopologues)

    @classmethod
    def from_dict(cls, isotopologues, temperature, scaling, imag_threshold,
                  reference_isotopologue="none",
                  mass_override_isotopologue="default", filename="<dict>"):
        """Build a Config programmatically, without a .config file.

        ``isotopologues`` maps each isotopologue name to a list of
        ``(gs_atom, ts_atom, replacement)`` tuples (1-indexed atoms).
        ``replacement`` is either a standard isotope label from weights.dat
        (e.g. ``"13C"``) or a number giving an unusual mass directly (e.g.
        ``5000.0``).
        """
        self = cls.__new__(cls)
        config = {i: None for i in EXPECTED_FIELDS}
        config["filename"] = filename
        config["temperature"] = temperature
        config["scaling"] = scaling
        config["imag_threshold"] = imag_threshold
        config["reference_isotopologue"] = reference_isotopologue
        config["mass_override_isotopologue"] = mass_override_isotopologue

        def _coerce(r):
            # a numeric string ("5000.0") is treated as a mass, matching the
            # .config file parser; a label ("13C") is left as-is
            if isinstance(r, str):
                try:
                    return float(r)
                except ValueError:
                    return r
            return r

        normalized = OrderedDict()
        for name, rules in isotopologues.items():
            if str(name) == "default":
                raise ValueError("name default is reserved.")
            normalized_rules = []
            for (from_atom, to_atom, replacement) in rules:
                if from_atom < 1 or to_atom < 1:
                    raise ValueError("check atom numbers for isotopologue %s" % name)
                replacement = _coerce(replacement)
                replacement_mass(replacement)   # validates label or numeric mass
                normalized_rules.append((int(from_atom), int(to_atom), replacement))
            normalized[str(name)] = normalized_rules

        self._finalize(config, normalized)
        return self

    def _finalize(self, config, isotopologues):
        # ensure we have all the fields we are supposed to
        for k,v in config.items():
            if k == "frequency_threshold" and v is not None:
                logger.warning("frequency_threshold is now deprecated and will be ignored.")
            if v is None:
                if k == "frequency_threshold":
                    config["frequency_threshold"] = 0.0
                else:
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

        config["reference_isotopologue"] = str(config["reference_isotopologue"])
        config["mass_override_isotopologue"] = str(config["mass_override_isotopologue"])

        config["frequency_threshold"] = float(config["frequency_threshold"])

        config["imag_threshold"] = float(config["imag_threshold"])
        if config["imag_threshold"] > 300.0:
            raise ValueError("imag threshold is too high")

        # store all the information in the object dictionary
        config["isotopologues"] = isotopologues
        self.__dict__ = config

    # checks if this config file is compatible with a pair of ground and transition state systems
    def check(self, gs, ts, verbose=False):
        # check that the isotopic replacements make sense
        isotopologues = self.isotopologues
        for i,isotopologue in isotopologues.items():
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
                # a custom numeric mass carries no element identity, so the
                # element-match check only applies to standard isotope labels
                if replacement not in REPLACEMENTS_Z:
                    continue
                replacementZ = REPLACEMENTS_Z[replacement]
                if from_atomZ != replacementZ:
                    raise ValueError("gs atomic number and replacement atomic number do not match for {0}\n".format(replacement_string))
                if to_atomZ != replacementZ:
                    raise ValueError("ts atomic number and replacement atomic number do not match for {0}\n".format(replacement_string))
                # unreachable: if we get here both gs and ts atomic numbers
                # already equal replacementZ, so they equal each other
                if (from_atomZ != to_atomZ):  # pragma: no cover
                    raise ValueError("gs and ts atomic number do not match for %s\n" % replacement_string)

        if verbose:
            logger.info("Config file %s makes sense with gs file %s and ts file %s.",
                        self.filename, gs.filename, ts.filename)

    # convert to human-readable format
    def __str__(self):
        to_string = "Config file: %s\nTemperature: %.1f K\nScaling: %.3f\nReference Isotopologue: %s\nImag threshold (cm-1): %d\n" % \
                    (self.filename, self.temperature, self.scaling, self.reference_isotopologue, self.imag_threshold)
        if self.frequency_threshold != 0:
            to_string += "Frequency threshold (cm-1): %d\n" % self.frequency_threshold

        keys = list(self.isotopologues.keys())
        if self.reference_isotopologue != "default" and self.reference_isotopologue != "none":
            try:
                keys.remove(self.reference_isotopologue)
            except ValueError:
                raise ValueError("Could not find the reference isotopologue: %s"
                                 % self.reference_isotopologue)
        keys.sort()

        if self.reference_isotopologue == "default" or self.reference_isotopologue == "none":
            to_string += "   No reference isotopologue.\n"
        else:
            keys = [self.reference_isotopologue] + keys
        for i in keys:
            isotopologue = self.isotopologues[i]

            for j in range(len(isotopologue)):
                to_string += "   Isotopologue {0: >10s}, replacement {1: >2d}: replace gs atom {2: ^3d} and ts atom {3: ^3d} with {4!s: >3}\n".format(i, j+1, isotopologue[j][0], isotopologue[j][1], isotopologue[j][2])

        return to_string[:-1]
