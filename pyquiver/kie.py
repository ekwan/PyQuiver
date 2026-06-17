# This calculates KIEs based on the Bigeleisen-Mayer equation.
import logging
import numpy as np
from . import quiver
from .config import Config
from .constants import DEFAULT_MASSES
from .results import Results, KIEResult, EIEResult
from collections import OrderedDict

# load constants
from .constants import PHYSICAL_CONSTANTS, replacement_mass

logger = logging.getLogger("pyquiver")
h  = PHYSICAL_CONSTANTS["h"]  # in J . s
c  = PHYSICAL_CONSTANTS["c"]  # in cm . s
kB = PHYSICAL_CONSTANTS["kB"] # in J/K

# a hydrogen substitution is treated as a *primary* H/D KIE when the substituted
# atom carries at least this fraction of the (mass-weighted) reaction mode, i.e.
# it is part of the reaction coordinate; secondary hydrogens carry far less
PRIMARY_HYDROGEN_MODE_FRACTION = 0.1

class KIE_Calculation(object):
    def __init__(self, config, gs, ts, style="gaussian", n_jobs=1):
        # n_jobs controls optional parallelism over isotopologues (default
        # serial). The heavy step is np.linalg.eigvalsh, which releases the
        # GIL, so threads give real speedup for large systems / many
        # isotopologues. n_jobs < 0 lets the pool pick a default worker count.
        self.n_jobs = n_jobs

        # accept either file paths (str) or already-built objects
        if isinstance(config, str):
            self.config = Config(config)
        elif isinstance(config, Config):
            self.config = config
        else:
            raise TypeError("config argument must be either a filepath or Config object.")

        if isinstance(gs, str):
            self.gs_system = quiver.System(gs, style=style)
        elif isinstance(gs, quiver.System):
            self.gs_system = gs
        else:
            raise TypeError("gs argument must be either a filepath or quiver.System object.")

        if isinstance(ts, str):
            self.ts_system = quiver.System(ts, style=style)
        elif isinstance(ts, quiver.System):
            self.ts_system = ts
        else:
            raise TypeError("ts argument must be either a filepath or quiver.System object.")

        # set the eie_flag to the recognized uninitialized value (used for checking if there are inconsistent calculation types)
        self.eie_flag = -1

        logger.debug("%s", self.config)

        KIES = OrderedDict()

        if self.config.frequency_threshold:
            logger.warning("config file uses the frequency_threshold parameter. This has been deprecated and low frequencies are dropped by linearity detection.")

        iso_pairs = list(self.make_isotopologues())

        # pre-compute the reference ("default") isotopologue's frequencies once,
        # serially: every pair shares those objects, so warming the cache here
        # both preserves the single-diagonalization optimization and avoids a
        # cache race when the per-isotopologue work runs on multiple threads
        if iso_pairs:
            (default_gs, _), (default_ts, _) = iso_pairs[0]
            default_gs.calculate_frequencies(self.config.imag_threshold,
                                             scaling=self.config.scaling)
            default_ts.calculate_frequencies(self.config.imag_threshold,
                                             scaling=self.config.scaling)

        def _build(pair):
            gs_tuple, ts_tuple = pair
            name = gs_tuple[1].name
            return name, KIE(name, gs_tuple, ts_tuple, self.config.temperature,
                             self.config.scaling, self.config.imag_threshold)

        if self.n_jobs == 1 or len(iso_pairs) <= 1:
            built = [_build(p) for p in iso_pairs]
        else:
            from concurrent.futures import ThreadPoolExecutor
            max_workers = self.n_jobs if self.n_jobs > 0 else None
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                built = list(pool.map(_build, iso_pairs))  # preserves order

        for name, kie in built:
            KIES[name] = kie
            self._warn_if_primary_hydrogen(name, kie)


        for name,k in KIES.items():
            if name != self.config.reference_isotopologue:
                if self.eie_flag == -1:
                    eie_flag_iso = name
                    self.eie_flag = k.eie_flag
                else:
                    # defensive: all isotopologues share the same TS, so they
                    # are all KIE or all EIE; a mismatch shouldn't occur
                    if self.eie_flag != k.eie_flag:  # pragma: no cover
                        if self.eie_flag == 1:
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
        keys = list(self.KIES.keys())
        
        # don't report the reference isotopologue (the mass-override one is
        # already excluded from KIES upstream)
        if self.config.reference_isotopologue in keys:
            keys.remove(self.config.reference_isotopologue)

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

        return (title_row, row, self.eie_flag)

    @property
    def results(self):
        """Structured, tabular view of the results (see pyquiver.results)."""
        return Results(self)

    def to_dict(self):
        """Map isotopologue name -> result value(s)."""
        return self.results.to_dict()

    def to_csv(self, path=None):
        """Render the results as CSV text (and optionally write to ``path``)."""
        return self.results.to_csv(path)

    def to_dataframe(self):
        """Return the results as a pandas DataFrame (optional pandas extra)."""
        return self.results.to_dataframe()

    def skodje_truhlar(self, reactant_energy, product_energy, ts_energy,
                       unit="hartree"):
        """Skodje-Truhlar tunnelling-corrected KIEs, one per isotopologue.

        Supply the single-point electronic energies of the reactant, product,
        and transition state (in hartree by default, or joules with
        ``unit='J'``); the barrier is V_ts - max(V_reactant, V_product). Each
        isotopologue's transition-state imaginary frequencies are taken from the
        calculation. Non-reference isotopologues are referenced exactly as in
        the printed table; the reference (if any) is returned as its absolute
        corrected KIE. Returns ``{name: corrected_kie}``.
        """
        from .tunneling import skodje_truhlar as _st_ratio

        if self.eie_flag != 0:
            raise ValueError("Skodje-Truhlar correction applies to KIE (not "
                             "EIE) calculations")
        if unit == "hartree":
            scale = PHYSICAL_CONSTANTS["Eh"]
        elif unit in ("J", "joule", "joules"):
            scale = 1.0
        else:
            raise ValueError("unit must be 'hartree' or 'J'")
        barrier = (ts_energy - max(reactant_energy, product_energy)) * scale

        def _imag(kie):
            light = kie.ts_tuple[0].calculate_frequencies(
                self.config.imag_threshold, scaling=self.config.scaling)[1]
            heavy = kie.ts_tuple[1].calculate_frequencies(
                self.config.imag_threshold, scaling=self.config.scaling)[1]
            return min(heavy), min(light)

        ref = self.config.reference_isotopologue
        has_ref = ref not in ("default", "none")
        if has_ref:
            h_heavy, h_light = _imag(self.KIES[ref])
            ref_ratio = _st_ratio(h_heavy, h_light, self.config.temperature, barrier)
        else:
            ref_ratio = 1.0

        corrected = OrderedDict()
        for name, kie in self.KIES.items():
            heavy, light = _imag(kie)
            ratio = _st_ratio(heavy, light, self.config.temperature, barrier)
            if name == ref and has_ref:
                corrected[name] = float(kie.value[0]) * ratio          # absolute
            else:
                corrected[name] = float(kie.value[0]) * ratio / ref_ratio
        return corrected

    def _warn_if_primary_hydrogen(self, name, kie):
        """Warn when an isotopologue is a primary hydrogen KIE, where the
        harmonic tunnelling corrections are unreliable. Treated as primary
        hydrogen when a substitution replaces a hydrogen (Z=1) that carries a
        large fraction of the transition-state reaction mode (i.e. the H is in
        the reaction coordinate)."""
        if kie.eie_flag != 0:   # no reaction mode for an EIE
            return
        znums = self.ts_system.atomic_numbers
        default_masses, sub_masses = kie.ts_tuple[0].masses, kie.ts_tuple[1].masses
        substituted_hydrogens = [i for i, (a, b) in enumerate(zip(default_masses, sub_masses))
                                 if a != b and znums[i] == 1]
        if not substituted_hydrogens:
            return
        composition = kie.ts_tuple[0].reaction_mode_composition()
        if any(composition[i] > PRIMARY_HYDROGEN_MODE_FRACTION for i in substituted_hydrogens):
            logger.warning("isotopologue %s looks like a primary hydrogen KIE "
                           "(the substituted H participates strongly in the "
                           "reaction mode); the harmonic Wigner and Bell "
                           "tunnelling corrections are unreliable here", name)

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
        for i in range(system.number_of_atoms):
            if not system.atomic_numbers[i] in DEFAULT_MASSES:
                raise ValueError("Default mass not available for atomic number %d at atom number %d in %s!" % (system.atomic_numbers[i], i+1, system.filename))
            masses.append(DEFAULT_MASSES[system.atomic_numbers[i]])
        return np.array(masses)

    def apply_mass_rules(self, prev_masses, rules):
        masses = list(prev_masses)
        for k,v in rules.items():
            masses[k] = v
        return masses

    def compile_mass_rules(self, iso_rules):
        gs_rules = {}
        ts_rules = {}
        for replacement in iso_rules:
            gs_atom_number, ts_atom_number, replacement_label = replacement
            mass = replacement_mass(replacement_label)
            gs_rules[gs_atom_number-1] = mass
            ts_rules[ts_atom_number-1] = mass
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

        for id_,iso in config.isotopologues.items():
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
        keys = list(self.KIES.keys())
        if self.config.reference_isotopologue in keys:
            keys.remove(self.config.reference_isotopologue)
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

        logger.debug("Calculating KIE for isotopologue %s.", name)
        self.value = self.calculate_kie()

    def calculate_kie(self):
        logger.debug("  Calculating Reduced Partition Function Ratio for Ground State.")
        rpfr_gs, gs_imag_ratios, gs_heavy_freqs, gs_light_freqs = calculate_rpfr(self.gs_tuple, self.imag_threshold, self.scaling, self.temperature)
        logger.debug("    rpfr_gs: %s", np.prod(rpfr_gs))
        logger.debug("  Calculating Reduced Partition Function Ratio for Transition State.")

        rpfr_ts, ts_imag_ratios, ts_heavy_freqs, ts_light_freqs = calculate_rpfr(self.ts_tuple, self.imag_threshold, self.scaling, self.temperature)
        logger.debug("    rpfr_ts: %s", np.prod(rpfr_ts))

        # eie_flag is -1 here (set once per KIE in __init__), so the else
        # branches below are defensive only
        if ts_imag_ratios is not None:
            if self.eie_flag == -1:
                self.eie_flag = 0
            else:  # pragma: no cover
                raise ValueError("quiver attempted to run a KIE calculation after an EIE calculation. Check the frequency threshold.")

            kies = ts_imag_ratios * rpfr_gs/rpfr_ts
            return kies
        else:
            if self.eie_flag == -1:
                self.eie_flag = 1
            else:  # pragma: no cover
                raise ValueError("quiver attempted to run a KIE calculation after an EIE calculation. Check the frequency threshold.")

            eie = rpfr_gs/rpfr_ts
            return eie

    def apply_reference(self, reference_kie):
        self.value /= reference_kie.value
        return self.value

    @property
    def result(self):
        """Named view of this isotopologue's result (KIEResult or EIEResult)."""
        if self.eie_flag == 1:
            return EIEResult(self.name, float(self.value))
        return KIEResult(self.name, float(self.value[0]),
                         float(self.value[1]), float(self.value[2]))

    def __str__(self):
        if self.value is not None:
            if self.eie_flag == 1:
                return "Isotopologue {1: >10s} {0: >33s} {2: ^12.4f} ".format("", self.name, self.value)
            else:
                return "Isotopologue {1: >10s} {0: >33s} {2: ^12.4f} {3: ^14.4f} {4: ^17.4f}".format("", self.name, self.value[0], self.value[1], self.value[2])
        else:  # pragma: no cover - value is always set in __init__
            return "KIE Object for isotopomer {0}. No value has been calculated yet.".format(self.name)

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
    # vectorized over all modes at once; u() is elementwise so numpy arrays
    # flow straight through (Wolfsberg eqn 4.79, one row per frequency)
    freqs_light = np.asarray(freqs_light, dtype=float)
    freqs_heavy = np.asarray(freqs_heavy, dtype=float)

    product_factor = freqs_heavy / freqs_light
    u_light = u(freqs_light, temperature)
    u_heavy = u(freqs_heavy, temperature)
    excitation_factor = (1.0 - np.exp(-u_light)) / (1.0 - np.exp(-u_heavy))
    ZPE_factor = np.exp(0.5 * (u_light - u_heavy))

    components = np.column_stack([product_factor, excitation_factor, ZPE_factor])

    if logger.isEnabledFor(logging.DEBUG):
        overall = product_factor * excitation_factor * ZPE_factor
        for i in range(len(freqs_light)):
            logger.debug("MODE %3d    LIGHT: %9.3f cm-1    HEAVY: %9.3f cm-1    FRQ RATIO: %9.5f    ZPE FACTOR: %9.5f    CONTRB TO RIPF: %9.5f",
                         i+1, freqs_light[i], freqs_heavy[i], product_factor[i], ZPE_factor[i], overall[i])

    return components

# tup is a tuple of a form (light_isotopologue, heavy_isotopologue)
def calculate_rpfr(tup, imag_threshold, scaling, temperature):
    # calculate_frequencies gives tuples of the form (small_freqs, imaginary_freqs, freqs)
    light_small_freqs, light_imag_freqs, light_freqs, light_num_small = tup[0].calculate_frequencies(imag_threshold, scaling=scaling)
    heavy_small_freqs, heavy_imag_freqs, heavy_freqs, heavy_num_small = tup[1].calculate_frequencies(imag_threshold, scaling=scaling)

    if len(heavy_freqs) != len(light_freqs):  # pragma: no cover - defensive
        raise ValueError("mismatch in the number of frequencies between isotopomers!")
    if len(light_imag_freqs) != len(heavy_imag_freqs):
        logger.warning("mismatch in the number of imaginary frequencies between isotopomers, ignoring imaginary mode")
        light_imag_freqs = []
        heavy_imag_freqs = []

    if logger.isEnabledFor(logging.DEBUG):
        def _fmt(values, spec="%.3f"):
            return "  ".join(spec % v for v in values) if len(values) else "none"
        logger.debug("light imaginary frequencies: %s", _fmt(light_imag_freqs))
        logger.debug("light small frequencies: %s", _fmt(light_small_freqs, "%.1f"))
        logger.debug("heavy imaginary frequencies: %s", _fmt(heavy_imag_freqs))
        logger.debug("heavy small frequencies: %s", _fmt(heavy_small_freqs, "%.1f"))

    raw_imag_ratio = None
    imag_ratios = None
    try:
        raw_imag_ratio = light_imag_freqs[0]/heavy_imag_freqs[0]
    except IndexError:
        pass

    if raw_imag_ratio is not None:
        wigner_imag_ratio = raw_imag_ratio * wigner(heavy_imag_freqs[0], light_imag_freqs[0], temperature)
        bell_imag_ratio = raw_imag_ratio * bell(heavy_imag_freqs[0], light_imag_freqs[0], temperature)
        imag_ratios = np.array([raw_imag_ratio, wigner_imag_ratio, bell_imag_ratio])

    partition_factors = partition_components(heavy_freqs, light_freqs, temperature)

    if logger.isEnabledFor(logging.DEBUG):
        factors = np.prod(partition_factors, axis=0)
        logger.debug("Product Factor: %s  Excitation Factor: %s  ZPE Factor: %s",
                     factors[0], factors[1], factors[2])

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
