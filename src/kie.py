# This calculates KIEs based on the Bigeleisen-Mayer equation.
import numpy as np

# load constants
from contants import PHYSICAL_CONSTANTS
h  = PHYSICAL_CONSTANTS["h"]  # in J . s
c  = PHYSICAL_CONSTANTS["c"]  # in cm . s
kB = PHYSICAL_CONSTANTS["kB"] # in J/K

# calculates the uncorrected EIE/KIE for a pair of isotopomers

def raw_KIE(gs_frequencies, ts_frequencies, temperature, frequency_threshold):
    # calculate the reduced isotopic partition function 
    pass

# calculates the reduced isotopic function ratio for a species (Wolfsberg eqn 4.79)
# assuming the symmetry ratio is 1/1
# uses the Teller-Redlich product rule
# returns 1 x n array of the partition function ratios, where n is the number of frequencies
# frequencies below frequency_threshold will be ignored and will not be included in the array
def partition_components(freqs_light, freqs_heavy, frequency_threshold, temperature):
    components = []
    for wavenumber_light, wavenumber_heavy in zip(freqs_light, freqs_heavy):
        if wavenumber_light < frequency_threshold:
            continue
        product_factor = wavenumber_heavy/wavenumber_light
        u_light = u(wavenumber_light, temperature)
        u_heavy = u(wavenumber_heavy, temperature)
        excitation_factor = (1.0-np.exp(-u_light))/(1.0-np.exp(-u_heavy))
        ZPE_factor = np.exp(0.5*(u_light-u_heavy))
        components.append(product_factor*excitation_factor*ZPE_factor)
    return np.array(components)

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
    else
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
    else
        correction_factor = 1.0
    return raw_KIEs * correction_factor


