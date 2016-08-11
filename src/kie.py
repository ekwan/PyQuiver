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
    def __init__(self, config_file, gs_file, ts_file, style="g09"):
        self.gs_system = System(gs_file, style=style)
        self.ts_system = System(ts_file, style=style)
        
        self.config = Config(config_file)

        print self.config.temperature
        KIES = {}
        for p in self.make_isotopologues():
            gs_tuple, ts_tuple = p
            name = gs_tuple[1].name
            print name
            kie = KIE(name, gs_tuple, ts_tuple, self.config.temperature, self.config)
            KIES[name] = kie
        
        for name,k in KIES.iteritems():
            if name != self.config.reference_isotopologue:
                if self.config.reference_isotopologue != "default":
                    k.apply_reference(KIES[self.config.reference_isotopologue])
                print k
        
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

        default_gs_masses = self.build_default_masses(self.gs_system)
        default_ts_masses = self.build_default_masses(self.ts_system)

        default_gs = Isotopologue("default", gs_system, default_gs_masses)
        default_ts = Isotopologue("default", ts_system, default_ts_masses)

        for id_,iso in config.isotopologues.iteritems():
            gs_rules, ts_rules = self.compile_mass_rules(iso)
            gs_masses = self.apply_mass_rules(default_gs_masses, gs_rules)
            ts_masses = self.apply_mass_rules(default_ts_masses, ts_rules)
            sub_gs = Isotopologue(id_, gs_system, gs_masses)
            sub_ts = Isotopologue(id_, ts_system, ts_masses)
            yield ((default_gs, sub_gs), (default_ts, sub_ts))
            #yield ((sub_gs, ref_gs), (sub_ts, ref_ts))
                
class KIE(object):
    # the constructor expects a tuple of the form yielded by make_isotopologue
    def __init__(self, name, gs_tuple, ts_tuple, temperature, config):
        # copy fields
        # the associated calculation object useful for pulling config fields etc.
        self.eie_flag = -1
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
        #print "Heavy frequencies:", freqs_heavy
        #print "Light frequencies:", freqs_light
        temperature = self.temperature
        components = []

        for wavenumber_light, wavenumber_heavy in zip(freqs_light, freqs_heavy):
            product_factor = wavenumber_heavy/wavenumber_light
            u_light = u(wavenumber_light, temperature)
            u_heavy = u(wavenumber_heavy, temperature)
            excitation_factor = (1.0-np.exp(-u_light))/(1.0-np.exp(-u_heavy))
            ZPE_factor = np.exp(0.5*(u_light-u_heavy))
            components.append([product_factor,excitation_factor,ZPE_factor])
        return np.array(components)

    # how should scaling factor be used?
    # tup is a tuple of a form 
    def calculate_rpfr(self, tup):
        # calculate_frequencies gives tuples of the form (small_freqs, imaginary_freqs, freqs)
        #print "Frequency threshold:", self.freq_threshold
        _, light_imag_freqs, light_freqs = tup[0].calculate_frequencies(self.freq_threshold, scaling=self.config.scaling)
        _, heavy_imag_freqs, heavy_freqs = tup[1].calculate_frequencies(self.freq_threshold, scaling=self.config.scaling)
        imag_ratio = None
        try:
            imag_ratio = light_imag_freqs[0]/heavy_imag_freqs[0]
            #print light_imag_freqs[0], heavy_imag_freqs[0]
            #print "IMAGINARY RATIO:", imag_ratio
        except:
            pass
            #print "No imaginary frequencies."
        '''
        print "Small frequencies:"
        print tup[0].frequencies[0]
        print tup[1].frequencies[0]
        print "Imaginary frequencies:"
        print tup[0].frequencies[1]
        print tup[1].frequencies[1]
        '''
        partition_factors = self.partition_components(heavy_freqs, light_freqs)
        return (np.prod(partition_factors), imag_ratio, heavy_freqs, light_freqs)

    def calculate_kie(self):

        rpfr_gs, gs_imag_ratio, gs_heavy_freqs, gs_light_freqs = self.calculate_rpfr(self.gs_tuple)
        #print "rpfr_gs:", np.prod(rpfr_gs)
        rpfr_ts, ts_imag_ratio, ts_heavy_freqs, ts_light_freqs = self.calculate_rpfr(self.ts_tuple)
        #print "rpfr_ts:", np.prod(rpfr_ts)

        if ts_imag_ratio:
            if self.eie_flag == -1:
                self.eie_flag = 0
            else:
                raise ValueError("quiver attempted to run a KIE calculation after an EIE calculation. Check the frequency threshold.")

            uncorrected_kie = ts_imag_ratio * rpfr_gs/rpfr_ts
            wigner_corrected = self.wigner(uncorrected_kie, ts_light_freqs, ts_heavy_freqs, self.config.temperature)
            bell_corrected = self.bell(uncorrected_kie, ts_light_freqs, ts_heavy_freqs)
            kie = np.array([uncorrected_kie, wigner_corrected, bell_corrected])
            return kie
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

    # calculates the Wigner tunnelling correction
    # multiplies the KIE by a factor of (1+u_H^2/24)/(1+u_D^2/24)
    # assumes the frequencies are sorted in ascending order
    def wigner(self, raw_KIE, ts_frequencies_light, ts_frequencies_heavy, temperature):
        # calculate correction factor
        wavenumber_light = ts_frequencies_light[0]
        wavenumber_heavy = ts_frequencies_heavy[0]
        if wavenumber_light < 0.0 and wavenumber_heavy < 0.0:
            u_H = u(wavenumber_light, temperature)
            u_D = u(wavenumber_heavy, temperature)
            correction_factor = (1.0 + u_H**2 / 24.0) / (1.0 + u_D**2 / 24.0)
        else:
            correction_factor = 1.0
        return raw_KIE * correction_factor

    # calculates the Bell infinite parabola tunnelling correction
    # multiplies the KIE by a factor of (u_H/u_D)*(sin(u_D/2)/sin(u_H/2))
    # assumes the frequencies are sorted in ascending order
    def bell(self, raw_KIEs, ts_frequencies_light, ts_frequencies_heavy):
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



    def __str__(self):
        return "isotopomer {0} {1}".format(self.name, self.value)


# utility function to calculate u terms
# u = hv/kT where h = Planck's constant, v = frequency, kB = Boltzmann's constant, T = temperature
# if using v in wavenumber, w
# w = v/c, with c in cm/s
# v = cw
# u = hcw/kT, with c in cm/s
def u(wavenumber, temperature):
    return h*c*wavenumber/(kB*temperature)

