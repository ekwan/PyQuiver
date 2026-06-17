"""Tunnelling corrections for transition-state imaginary modes.

This module gathers the three corrections PyQuiver supports:

* :func:`wigner` and :func:`bell` (re-exported from :mod:`pyquiver.kie`) need
  only the imaginary frequencies and the temperature.
* :func:`skodje_truhlar` additionally needs the reaction barrier height
  (J. Phys. Chem. 1981, 85, 624-626). It is more appropriate than the inverted
  parabola for very large imaginary frequencies.

All three take the imaginary frequencies as negative wavenumbers (cm^-1), the
convention PyQuiver uses internally, and return the *ratio* of the light to the
heavy correction (i.e. the factor a KIE is multiplied by).
"""

import numpy as np

from .constants import PHYSICAL_CONSTANTS
from .kie import wigner, bell  # re-export so all corrections share one namespace

h = PHYSICAL_CONSTANTS["h"]    # J s
c = PHYSICAL_CONSTANTS["c"]    # cm / s
kB = PHYSICAL_CONSTANTS["kB"]  # J / K
Eh = PHYSICAL_CONSTANTS["Eh"]  # J / hartree

__all__ = ["wigner", "bell", "skodje_truhlar", "skodje_truhlar_kappa"]


def skodje_truhlar_kappa(imag_wavenumber, temperature, barrier):
    """Skodje-Truhlar transmission coefficient for a single imaginary mode.

    ``barrier`` is the barrier height in joules (V = V_ts - max(V_reactant,
    V_product)). ``imag_wavenumber`` is the imaginary mode in cm^-1 (sign is
    ignored; its magnitude is the natural frequency).
    """
    imag_hz = np.abs(imag_wavenumber * c)
    alpha = 2.0 * np.pi / (h * imag_hz)
    beta = 1.0 / (kB * temperature)
    bpa = beta * np.pi / alpha

    if alpha > beta:
        return bpa / np.sin(bpa) - beta / (alpha - beta) * np.exp((beta - alpha) * barrier)
    elif alpha < beta:
        return beta / (beta - alpha) * (np.exp((beta - alpha) * barrier) - 1.0)
    else:  # pragma: no cover - alpha == beta is a measure-zero edge case
        return bpa / np.sin(bpa)


def skodje_truhlar(ts_imag_heavy, ts_imag_light, temperature, barrier):
    """Skodje-Truhlar tunnelling correction: kappa_light / kappa_heavy.

    ``ts_imag_heavy``/``ts_imag_light`` are the heavy/light transition-state
    imaginary frequencies (negative cm^-1); ``barrier`` is in joules.
    """
    if not (ts_imag_heavy < 0.0 and ts_imag_light < 0.0):
        raise ValueError("imaginary frequency passed to Skodje-Truhlar correction was real")
    light = skodje_truhlar_kappa(ts_imag_light, temperature, barrier)
    heavy = skodje_truhlar_kappa(ts_imag_heavy, temperature, barrier)
    return light / heavy
