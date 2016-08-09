# This calculates KIEs based on the Bigeleisen-Mayer equation.
import numpy as np
from quiver import System, Isotopologue
from config import Config
from constants import DEFAULT_MASSES

# load constants
from constants import PHYSICAL_CONSTANTS, REPLACEMENTS
h  = PHYSICAL_CONSTANTS["h"]  # in J . s
c  = PHYSICAL_CONSTANTS["c"]  # in cm . s
kB = PHYSICAL_CONSTANTS["kB"] # in J/K

class KIE_Calculation(object):
    def __init__(self, config_file, gs_file, ts_file):
        self.gs_system = System(gs_file)
        self.ts_system = System(ts_file)
        
        self.config = Config(config_file)

        for p in self.make_isotopologues():
            gs_tuple, ts_tuple = p
            kie = KIE(gs_tuple[0].name, gs_tuple, ts_tuple, self.config.temperature, self.config)
            print kie
            #print p[0][0].name
            #print "Sub masses:", p[0][0].masses
            #print p[0][1].name
            #print "Ref masses:", p[0][1].masses
        
    def build_reference_masses(self):
        config = self.config
        gs_system = self.gs_system
        ts_system = self.ts_system

        iso = self.config.isotopologues[config.reference_isotopologue]
        gs_rules, ts_rules = self.compile_mass_rules(iso)

        gs_masses = self.build_default_masses(gs_system)
        ts_masses = self.build_default_masses(ts_system)

        gs_masses = self.apply_mass_rules(gs_masses, gs_rules)
        ts_masses = self.apply_mass_rules(ts_masses, ts_rules)

        return gs_masses, ts_masses

    def build_default_masses(self, system):
        masses = []
        for i in xrange(system.number_of_atoms):
            masses.append(DEFAULT_MASSES[system.atomic_numbers[i]])
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

        ref_gs_masses, ref_ts_masses = self.build_reference_masses()

        ref_gs = Isotopologue(self.config.reference_isotopologue, gs_system, ref_gs_masses)
        ref_ts = Isotopologue(self.config.reference_isotopologue, ts_system, ref_ts_masses)

        for id_,iso in config.isotopologues.iteritems():
            if id_ == config.reference_isotopologue:
                pass
            else:
                gs_rules, ts_rules = self.compile_mass_rules(iso)
                gs_masses = self.apply_mass_rules(ref_gs_masses, gs_rules)
                ts_masses = self.apply_mass_rules(ref_ts_masses, ts_rules)
                sub_gs = Isotopologue(id_, gs_system, gs_masses)
                sub_ts = Isotopologue(id_, ts_system, ts_masses)
                
                yield ((sub_gs, ref_gs), (sub_ts, ref_ts))
                

class KIE(object):
    # the constructor expects a tuple of the form yielded by make_isotopologue
    def __init__(self, name, gs_tuple, ts_tuple, temperature, config):
        # copy fields
        # the associated calculation object useful for pulling config fields etc.
        self.name = name
        self.freq_threshold = config.frequency_threshold
        self.config = config
        self.gs_tuple, self.ts_tuple = gs_tuple, ts_tuple
        self.temperature = temperature

        self.value = self.calculate_kie()

    # calculates the reduced isotopic function ratio for a species (Wolfsberg eqn 4.79)
    # assuming the symmetry ratio is 1/1
    # uses the Teller-Redlich product rule
    # returns 1 x n array of the partition function ratios, where n is the number of frequencies
    # frequencies below frequency_threshold will be ignored and will not be included in the array
    def partition_components(self, freqs_heavy, freqs_light):
        print "Heavy frequencies:", freqs_heavy
        print "Light frequencies:", freqs_light
        temperature = self.temperature
        components = []
        for wavenumber_light, wavenumber_heavy in zip(freqs_light, freqs_heavy):
            product_factor = wavenumber_heavy/wavenumber_light
            u_light = u(wavenumber_light, temperature)
            u_heavy = u(wavenumber_heavy, temperature)
            excitation_factor = (1.0-np.exp(-u_light))/(1.0-np.exp(-u_heavy))
            ZPE_factor = np.exp(0.5*(u_light-u_heavy))
            components.append(product_factor*excitation_factor*ZPE_factor)
        return np.array(components)

    # how should scaling factor be used?
    def calculate_rpfr(self, tup):
        # calculate_frequencies gives tuples of the form (small_freqs, imaginary_freqs, freqs)
        light_freqs = tup[0].calculate_frequencies(self.freq_threshold)[2]
        heavy_freqs = tup[1].calculate_frequencies(self.freq_threshold)[2]
        
        return self.partition_components(heavy_freqs, light_freqs)

    def calculate_kie(self):

        rpfr_gs = self.calculate_rpfr(self.gs_tuple)
        rpfr_ts = self.calculate_rpfr(self.ts_tuple)

        kie = np.prod(rpfr_gs)/np.prod(rpfr_ts)
        return kie

    def __str__(self):
        return "isotopomer {0} {1}".format(self.name, self.value)

'''
# calculates the EIE or KIE for the specified systems
# returns a KIECalculation
def calculate_KIE(config, gs_system, ts_system):
    # get fields
    temperature = config.temperature
    scaling = config.scaling
    frequency_threshold = config.frequency_threshold
    reference_isotopologue = config.reference_isotopologue
    
    # generate isotopologues
    isotopologues = make_isotopologues(config, gs_system, ts_system)

    # calculate isotope effects
    for gs_isotopologue, ts_isotopologue, description in isotopologues:


    # construct result


# calculates the uncorrected EIE/KIE for a pair of isotopologues
# assumes the frequencies are sorted in ascending order
def uncorrected_KIE(gs_isotopologue_light, gs_isotopologue_heavy, ts_isotopologue_light, ts_isotopologue_heavy, \
                    temperature, frequency_threshold, frequency_method="mass weighted hessian"):
    # get the frequencies
    gs_freqs_light = gs_isotopologue

    if ( gs_frequencies[0] < 0.0 and ts_frequencies[0] < 0.0):
        # do a KIE calculation

        # add the tunnelling corrections

        # return result
    else:
        # do an EIE calculation
        partition_function_

        # return result
'''
# utility function to calculate u terms
# u = hv/kT where h = Planck's constant, v = frequency, kB = Boltzmann's constant, T = temperature
# if using v in wavenumber, w
# w = v/c, with c in cm/s
# v = cw
# u = hcw/kT, with c in cm/s
def u(wavenumber, temperature):
    return h*c*wavenumber/(kB*temperature)

# calculates the Wigner tunnelling correction
# multiplies the KIE by a factor of (1+u_H^2/24)/(1+u_D^2/24)
# assumes the frequencies are sorted in ascending order
def wigner(raw_KIEs, ts_frequencies_light, ts_frequencies_heavy, temperature):
    # calculate correction factor
    wavenumber_light = ts_frequencies_light[0]
    wavenumber_heavy = ts_frequencies_heavy[0]
    if wavenumber_light < 0.0 and wavenumber_heavy < 0.0:
        u_H = u(wavenumber_light, temperature)
        u_D = u(wavenumber_heavy, temperature)
        correction_factor = (1.0 + u_H**2 / 24.0) / (1.0 + u_D**2 / 24.0)
    else:
        correction_factor = 1.0
    return raw_KIEs * correction_factor

# calculates the Bell infinite parabola tunnelling correction
# multiplies the KIE by a factor of (u_H/u_D)*(sin(u_D/2)/sin(u_H/2))
# assumes the frequencies are sorted in ascending order
def bell(raw_KIEs, ts_frequencies_light, ts_frequencies_heavy):
    # calculate correction factor
    wavenumber_light = ts_frequencies_light[0]
    wavenumber_heavy = ts_frequencies_heavy[0]
    if wavenumber_light < 0.0 and wavenumber_heavy < 0.0:
        u_H = u(wavenumber_light, temperature)
        u_D = u(wavenumber_heavy, temperature)
        correction_factor = (u_H/u_D)*(np.sin(u_D/2.0)/np.sin(u_H/2.0))
    else:
        correction_factor = 1.0
    return raw_KIEs * correction_factor

