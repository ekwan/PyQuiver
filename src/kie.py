# This calculates KIEs based on the Bigeleisen-Mayer equation.
import numpy as np
from contants import PHYSICAL_CONSTANTS

# calculates the KIE for a pair of isotopomers
# returns a tuple of KIEs: (uncorrected, Wigner-corrected, infinite-parabola corrected) 
# for EIEs, we use Wolfsberg eqn 4.80
# for KIEs, we use Wolfsberg eqn 4.142
# Wigner correction: Wolfsberg eqn 4.137
def raw_KIE(gs_frequencies, ts_frequencies, temperature, frequency_threshold):
    pass

# calculates the reduced isotopic function ratio (Wolfsberg eqn 4.78)
# assuming the symmetry ratio is 1/1
# uses the Teller-Redlich rule
def partition(masses_unsubstituted, masses_substituted, frequencies, temperature, frequency_threshold):
    pass

# u = hv/kT
# if using v in wavenumber, w
# w = v/c, with c in cm/s
# v = cw
# u = hcw/kT, with c in cm/s
def u(frequency, temperature):
    pass

# calculates the Wigner tunnelling correction
def wigner():
    pass

# calculates the Bell infinite parabola tunnelling correction
def bell():
    pass
