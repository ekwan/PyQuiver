# This calculates KIEs based on the Bigeleisen-Mayer equation.
import numpy as np
#from quiver import System, Isotopologue, DEBUG
import quiver
import settings
from config import Config
from constants import DEFAULT_MASSES
from collections import OrderedDict

# load constants
from constants import PHYSICAL_CONSTANTS, REPLACEMENTS
h  = PHYSICAL_CONSTANTS["h"]  # in J . s
c  = PHYSICAL_CONSTANTS["c"]  # in cm . s
kB = PHYSICAL_CONSTANTS["kB"] # in J/K

class KIE_Calculation(object):
    def __init__(self, config, gs, ts, style="g09"):
        # check the types of config, gs, and ts parsing files if necessary and copying fields if not
        if type(config) is str:
            self.config = Config(config)
        elif type(config) is Config:
            self.config = Config
        else:
            raise TypeError("config argument must be either a filepath or Config object.")

        if type(gs) is str:
            self.gs_system = quiver.System(gs, style=style)
        elif type(gs) is quiver.System:
            self.gs_system = gs
        else:
            raise TypeError("gs argument must be either a filepath or quiver.System object.")

        if type(ts) is str:
            self.ts_system = quiver.System(ts, style=style)
        elif type(ts) is quiver.System:
            self.ts_system = ts
        else:
            raise TypeError("ts argument must be either a filepath or quiver.System object.")

        # set the eie_flag to the recognized uninitialized value (used for checking if there are inconsistent calculation types)
        self.eie_flag = -1

        if settings.DEBUG != 0:
            print self.config
            
        KIES = OrderedDict()

        if self.config.frequency_threshold:
            print "WARNING: config file uses the frequency_threshold parameter. This has been deprecated and low frequencies are dropped by linearity detection."

        for p in self.make_isotopologues():
            gs_tuple, ts_tuple = p
            name = gs_tuple[1].name
            kie = KIE(name, gs_tuple, ts_tuple, self.config.temperature, self.config.scaling, self.config.imag_threshold)
            KIES[name] = kie
        
        for name,k in KIES.iteritems():
            if name != self.config.reference_isotopologue:
                if self.eie_flag == -1:
                    eie_flag_iso = name
                    self.eie_flag = k.eie_flag
                else:
                    if self.eie_flag != k.eie_flag:
                        if eie_flag == 1:
                            raise ValueError("quiver attempted to run a KIE calculation (isotopomer {0}) after an EIE calculation (isotopomer {1}). Check the frequency threshold.".format(name, eie_flag_iso))
                        else:
                            raise ValueError("quiver attempted to run an EIE calculation (isotopomer {0}) after a KIE calculation (isotopomer {1}). Check the frequency threshold.".format(name, eie_flag_iso))

                if self.config.reference_isotopologue != "default" and self.config.reference_isotopologue != "none":
                    k.apply_reference(KIES[self.config.reference_isotopologue])

        self.KIES = KIES

    # retrieves KIEs for autoquiver output
    # if report_tunnelling = True, the first number will be the inverted parabola KIE
    # and the second number will be the tunnelling correction
    def get_row(self, report_tunnelling=False):
        title_row = ""
        row = ""
        keys = self.KIES.keys()
        
        # don't report the reference isotoplogue
        if self.config.reference_isotopologue != "default" and self.config.reference_isotopologue != "none":
            keys.remove(self.config.reference_isotopologue)

        if self.config.mass_override_isotopologue != "default" and self.config.reference_isotopologue != "none":
            keys.remove(self.config.mass_override_isotopologue)

        # for each isotopologue, add the KIEs
        for name in keys:
            if report_tunnelling:
                title_row += "%s,%s," % (name + "_uncorr", name + "_inf_para")
            else:
                title_row += "{0},".format(name)

            if self.eie_flag == 0:
                # this is a KIE calculation
                if report_tunnelling:
                    row += "%.4f,%.4f," % ( self.KIES[name].value[-1], self.KIES[name].value[-1]-self.KIES[name].value[-3] )
                else:
                    row += "{0:.4f},".format(self.KIES[name].value[-1])
            else:
                # this is an EIE calculation
                row += "{0:.4f},".format(self.KIES[name].value)
            #if settings.DEBUG >= 2:
            #    if self.eie_flag == 0:
            #        print "KIE calculation detected using {0}th tunneling correction".format(len(self.KIES[name].value))
            #        row += "{0:.3f},".format(self.KIES[name].value[-1])
            #    else:
            #        row += "{0:.3f},".format(self.KIES[name].value)

        return (title_row, row, self.eie_flag)

    def build_mass_override_masses(self):
        config = self.config
        gs_system = self.gs_system
        ts_system = self.ts_system

        gs_masses = self.build_default_masses(gs_system)
        ts_masses = self.build_default_masses(ts_system)

        if config.mass_override_isotopologue != "default":
            iso = self.config.isotopologues[config.mass_override_isotopologue]
            gs_rules, ts_rules = self.compile_mass_rules(iso)

            gs_masses = self.apply_mass_rules(gs_masses, gs_rules)
            ts_masses = self.apply_mass_rules(ts_masses, ts_rules)

        return gs_masses, ts_masses

    def build_default_masses(self, system):
        masses = []
        #print "---start---"
        #for k in DEFAULT_MASSES:
        #    print k, DEFAULT_MASSES[k]
        #print "---"
        for i in xrange(system.number_of_atoms):
            #print i, system.atomic_numbers[i], DEFAULT_MASSES[system.atomic_numbers[i]]
            if not system.atomic_numbers[i] in DEFAULT_MASSES:
                raise ValueError("Default mass not available for atomic number %d at atom number %d in %s!" % (system.atomic_numbers[i], i+1, system.filename))
            masses.append(DEFAULT_MASSES[system.atomic_numbers[i]])
        #print "--end--"
        return np.array(masses)

    def apply_mass_rules(self, prev_masses, rules):
        masses = list(prev_masses)
        for k,v in rules.iteritems():
            masses[k] = v
        return masses

    def compile_mass_rules(self, iso_rules):
        gs_rules = {}
        ts_rules = {}
        for replacement in iso_rules:
            gs_atom_number, ts_atom_number, replacement_label = replacement
            replacement_mass = REPLACEMENTS[replacement_label]
            gs_rules[gs_atom_number-1] = replacement_mass
            ts_rules[ts_atom_number-1] = replacement_mass
        return gs_rules, ts_rules

    # make the requested isotopic substitutions
    # yields tuples of tuples of the form ((gs_sub, gs_ref), (ts_sub, ts_ref))
    def make_isotopologues(self):
        config = self.config
        gs_system = self.gs_system
        ts_system = self.ts_system
        config.check(gs_system, ts_system)

        mass_override_gs_masses, mass_override_ts_masses = self.build_mass_override_masses()

        default_gs = quiver.Isotopologue("default", gs_system, mass_override_gs_masses)
        default_ts = quiver.Isotopologue("default", ts_system, mass_override_ts_masses)

        for id_,iso in config.isotopologues.iteritems():
            if id_ != config.mass_override_isotopologue:
                gs_rules, ts_rules = self.compile_mass_rules(iso)
                gs_masses = self.apply_mass_rules(mass_override_gs_masses, gs_rules)
                ts_masses = self.apply_mass_rules(mass_override_ts_masses, ts_rules)
                sub_gs = quiver.Isotopologue(id_, gs_system, gs_masses)
                sub_ts = quiver.Isotopologue(id_, ts_system, ts_masses)
                yield ((default_gs, sub_gs), (default_ts, sub_ts))
               
    def __str__(self):
        string = "\n=== PyQuiver Analysis ===\n"
        if self.eie_flag == 0:
            string += "Isotopologue                                              uncorrected      Wigner     inverted parabola\n"
            string += "                                                              KIE           KIE              KIE"
        else:
            string += "Isotopologue                                                  EIE"
        keys = self.KIES.keys()
        if self.config.reference_isotopologue != "default" and self.config.reference_isotopologue != "none":
            keys.remove(self.config.reference_isotopologue)
        if self.config.mass_override_isotopologue != "default":
            keys.remove(self.config.mass_override_isotopologue)
        #keys.sort()
        for name in keys:
            string += "\n" + str(self.KIES[name])

        if self.config.reference_isotopologue != "default" and self.config.reference_isotopologue != "none":
            string += "\n\nKIEs referenced to isotopologue {0}, whose absolute KIEs are:".format(self.config.reference_isotopologue)
            string += "\n" + str(self.KIES[self.config.reference_isotopologue])
        else:
            string += "\n\nAbsolute KIEs are given."

        return string
            
class KIE(object):
    # the constructor expects a tuple of the form yielded by make_isotopologue
    def __init__(self, name, gs_tuple, ts_tuple, temperature, scaling, imag_threshold):
        # copy fields
        # the associated calculation object useful for pulling config fields etc.
        self.eie_flag = -1
        self.name = name
        self.imag_threshold = imag_threshold
        self.gs_tuple, self.ts_tuple = gs_tuple, ts_tuple
        self.temperature = temperature
        self.scaling = scaling
        
        if settings.DEBUG >= 2:
            print "Calculating KIE for isotopologue {0}.".format(name)
        self.value = self.calculate_kie()

    def calculate_kie(self):
        if settings.DEBUG >= 2:
            print "  Calculating Reduced Partition Function Ratio for Ground State."        
        rpfr_gs, gs_imag_ratios, gs_heavy_freqs, gs_light_freqs = calculate_rpfr(self.gs_tuple, self.imag_threshold, self.scaling, self.temperature)
        if settings.DEBUG >= 2:
            print "    rpfr_gs:", np.prod(rpfr_gs)
        if settings.DEBUG >= 2:
            print "  Calculating Reduced Partition Function Ratio for Transition State."

        rpfr_ts, ts_imag_ratios, ts_heavy_freqs, ts_light_freqs = calculate_rpfr(self.ts_tuple, self.imag_threshold, self.scaling, self.temperature)
        if settings.DEBUG >= 2:
            print "    rpfr_ts:", np.prod(rpfr_ts)

        if ts_imag_ratios is not None:
            if self.eie_flag == -1:
                self.eie_flag = 0
            else:
                raise ValueError("quiver attempted to run a KIE calculation after an EIE calculation. Check the frequency threshold.")

            kies = ts_imag_ratios * rpfr_gs/rpfr_ts
            return kies
        else:
            if self.eie_flag == -1:
                self.eie_flag = 1
            else:
                raise ValueError("quiver attempted to run a KIE calculation after an EIE calculation. Check the frequency threshold.")
            
            eie = rpfr_gs/rpfr_ts
            return eie

    def apply_reference(self, reference_kie):
        self.value /= reference_kie.value
        return self.value

    def __str__(self):
        if self.value is not None:
            if self.eie_flag == 1:
                return "Isotopologue {1: >10s} {0: >33s} {2: ^12.4f} ".format("", self.name, self.value)
            else:
                return "Isotopologue {1: >10s} {0: >33s} {2: ^12.4f} {3: ^14.4f} {4: ^17.4f}".format("", self.name, self.value[0], self.value[1], self.value[2])
        else:
            "KIE Object for isotopomer {0}. No value has been calculated yet.".format(self.name)

# utility function to calculate u terms
# u = hv/kT where h = Planck's constant, v = frequency, kB = Boltzmann's constant, T = temperature
# if using v in wavenumber, w
# w = v/c, with c in cm/s
# v = cw
# u = hcw/kT, with c in cm/s
def u(wavenumber, temperature):
    return h*c*wavenumber/(kB*temperature)

# calculates the reduced isotopic function ratio for a species (Wolfsberg eqn 4.79)
# assuming the symmetry ratio is 1/1
# uses the Teller-Redlich product rule
# returns 1 x n array of the partition function ratios, where n is the number of frequencies
# frequencies below frequency_threshold will be ignored and will not be included in the array
def partition_components(freqs_heavy, freqs_light, temperature):
    components = []
    i = 0
    for wavenumber_light, wavenumber_heavy in zip(freqs_light, freqs_heavy):
        i += 1
        product_factor = wavenumber_heavy/wavenumber_light
        u_light = u(wavenumber_light, temperature)
        u_heavy = u(wavenumber_heavy, temperature)
        excitation_factor = (1.0-np.exp(-u_light))/(1.0-np.exp(-u_heavy))
        ZPE_factor = np.exp(0.5*(u_light-u_heavy))
        components.append([product_factor,excitation_factor,ZPE_factor])
        if settings.DEBUG >= 3:
            overall_factor = product_factor * excitation_factor * ZPE_factor
            print "MODE %3d    LIGHT: %9.3f cm-1    HEAVY: %9.3f cm-1    FRQ RATIO: %9.5f    ZPE FACTOR: %9.5f    CONTRB TO RIPF: %9.5f" % (i, wavenumber_light, wavenumber_heavy, product_factor, ZPE_factor, overall_factor)
            #if overall_factor < 0.99 or overall_factor > 1.01:
            #    print " *"
            #else:
            #    print
    return np.array(components)

# tup is a tuple of a form (light_isotopologue, heavy_isotopologue)
def calculate_rpfr(tup, imag_threshold, scaling, temperature):
    # calculate_frequencies gives tuples of the form (small_freqs, imaginary_freqs, freqs)
    light_small_freqs, light_imag_freqs, light_freqs, light_num_small = tup[0].calculate_frequencies(imag_threshold, scaling=scaling)
    heavy_small_freqs, heavy_imag_freqs, heavy_freqs, heavy_num_small = tup[1].calculate_frequencies(imag_threshold, scaling=scaling)
    
    if len(heavy_freqs) != len(light_freqs):
        raise ValueError("mismatch in the number of frequencies between isotopomers!")
    if len(light_imag_freqs) != len(heavy_imag_freqs):
        print "WARNING: mismatch in the number of imaginary frequencies between isotopomers, ignoring imaginary mode"
        light_imag_freqs = []
        heavy_imag_freqs = []

    if settings.DEBUG > 2:
        print "light imaginary frequencies: ",
        if len(light_imag_freqs) == 0:
            print "none",
        for i in light_imag_freqs:
            print "%.3f  " % i,
        print
        print "light small frequencies: ",
        for i in light_small_freqs:
            print "%.1f  " % i,
        print
        print "heavy imaginary frequencies: ",
        if len(heavy_imag_freqs) == 0:
            print "none",
        for i in heavy_imag_freqs:
            print "%.3f  " % i,
        print
        print "heavy small frequencies: ",
        for i in heavy_small_freqs:
            print "%.1f  " % i,
        print
   
    raw_imag_ratio = None
    imag_ratios = None
    try:
        raw_imag_ratio = light_imag_freqs[0]/heavy_imag_freqs[0]
    except IndexError:
        pass

    if raw_imag_ratio:
        wigner_imag_ratio = raw_imag_ratio * wigner(heavy_imag_freqs[0], light_imag_freqs[0], temperature)
        bell_imag_ratio = raw_imag_ratio * bell(heavy_imag_freqs[0], light_imag_freqs[0], temperature)
        imag_ratios = np.array([raw_imag_ratio, wigner_imag_ratio, bell_imag_ratio])

    partition_factors = partition_components(heavy_freqs, light_freqs, temperature)

    if settings.DEBUG >= 2:
        factors = np.prod(partition_factors, axis=0)
        print "{3: ^8}Product Factor: {0}\n{3: ^8}Excitation Factor: {1}\n{3: ^8}ZPE Factor: {2}".format(factors[0], factors[1], factors[2], "")

    return (np.prod(partition_factors), imag_ratios, np.array(heavy_freqs), np.array(light_freqs))

# calculates the Wigner tunnelling correction
# multiplies the KIE by a factor of (1+u_H^2/24)/(1+u_D^2/24)
# assumes the frequencies are sorted in ascending order
def wigner(ts_imag_heavy, ts_imag_light, temperature):
    # calculate correction factor
    if ts_imag_heavy < 0.0 and ts_imag_light < 0.0:
        u_H = u(ts_imag_light, temperature)
        u_D = u(ts_imag_heavy, temperature)
        correction_factor = (1.0 + u_H**2 / 24.0) / (1.0 + u_D**2 / 24.0)
    else:
        raise ValueError("imaginary frequency passed to Wigner correction was real")

    return correction_factor

# calculates the Bell inverted parabola tunneling correction
# multiplies the KIE by a factor of (u_H/u_D)*(sin(u_D/2)/sin(u_H/2))
# assumes the frequencies are sorted in ascending order
def bell(ts_imag_heavy, ts_imag_light, temperature):
    # calculate correction factor
    if ts_imag_heavy < 0.0 and ts_imag_light < 0.0:
        u_H = u(ts_imag_light, temperature)
        u_D = u(ts_imag_heavy, temperature)
        correction_factor = (u_H/u_D)*(np.sin(u_D/2.0)/np.sin(u_H/2.0))
    else:
        raise ValueError("imaginary frequency passed to Bell correction was real")
    return correction_factor
