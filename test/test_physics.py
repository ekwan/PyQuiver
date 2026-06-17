"""Unit tests for the pure physics functions in ``kie.py``:
the reduced vibrational term ``u``, the Wigner and Bell tunnelling
corrections, and the Bigeleisen-Mayer partition components.
"""

import numpy as np
import pytest

from pyquiver.kie import u, wigner, bell, partition_components, h, c, kB


# --- u = h c w / (kB T) ------------------------------------------------------

def test_u_matches_definition():
    assert u(1000.0, 300.0) == pytest.approx(h * c * 1000.0 / (kB * 300.0))


def test_u_zero_frequency():
    assert u(0.0, 300.0) == 0.0


def test_u_linear_in_wavenumber():
    assert u(2000.0, 300.0) == pytest.approx(2.0 * u(1000.0, 300.0))


def test_u_inverse_in_temperature():
    assert u(1000.0, 600.0) == pytest.approx(0.5 * u(1000.0, 300.0))


def test_u_representative_magnitude():
    # a ~1000 cm-1 mode at room temperature has u of order a few
    assert u(1000.0, 300.0) == pytest.approx(4.796, abs=1e-2)


# --- Wigner tunnelling correction -------------------------------------------

def test_wigner_unity_for_equal_frequencies():
    assert wigner(-500.0, -500.0, 300.0) == pytest.approx(1.0)


def test_wigner_lighter_isotope_tunnels_more():
    # the light isotope has the more negative imaginary frequency, so the ratio
    # (light correction / heavy correction) should exceed 1
    assert wigner(-400.0, -500.0, 300.0) > 1.0


def test_wigner_value_against_formula():
    u_H = u(-500.0, 300.0)
    u_D = u(-400.0, 300.0)
    expected = (1.0 + u_H ** 2 / 24.0) / (1.0 + u_D ** 2 / 24.0)
    assert wigner(-400.0, -500.0, 300.0) == pytest.approx(expected)


def test_wigner_rejects_real_frequencies():
    # tunnelling only applies to imaginary (negative) modes
    with pytest.raises(ValueError):
        wigner(500.0, 500.0, 300.0)


# --- Bell inverted-parabola correction --------------------------------------

def test_bell_unity_for_equal_frequencies():
    assert bell(-500.0, -500.0, 300.0) == pytest.approx(1.0)


def test_bell_value_against_formula():
    u_H = u(-500.0, 300.0)
    u_D = u(-400.0, 300.0)
    expected = (u_H / u_D) * (np.sin(u_D / 2.0) / np.sin(u_H / 2.0))
    assert bell(-400.0, -500.0, 300.0) == pytest.approx(expected)


def test_bell_rejects_real_frequencies():
    with pytest.raises(ValueError):
        bell(500.0, 500.0, 300.0)


# --- Bigeleisen-Mayer partition components ----------------------------------

def test_partition_components_shape():
    comps = partition_components([1100.0, 1200.0], [1000.0, 1100.0], 300.0)
    assert comps.shape == (2, 3)


def test_partition_components_identical_isotopologues_give_unity():
    freqs = [1000.0, 1500.0, 2000.0]
    comps = partition_components(freqs, freqs, 300.0)
    assert np.allclose(comps, 1.0)
    assert float(np.prod(comps)) == pytest.approx(1.0)


def test_partition_components_known_single_mode():
    heavy, light, T = 1000.0, 1100.0, 300.0
    u_l, u_h = u(light, T), u(heavy, T)
    expected = [heavy / light,
                (1.0 - np.exp(-u_l)) / (1.0 - np.exp(-u_h)),
                np.exp(0.5 * (u_l - u_h))]
    comps = partition_components([heavy], [light], T)
    assert comps[0] == pytest.approx(expected)
