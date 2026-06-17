"""Tests for mass-weighting, frequency calculation, and the caching
behaviour of the reference isotopologue.

The caching test is the regression guard for the performance discussion:
when several isotopologues share one config, the reference ("default")
isotopologue must be diagonalized only once, not once per substitution.
"""

from unittest import mock

import numpy as np
import pytest

from pyquiver import quiver
from pyquiver.kie import KIE_Calculation


class FakeSystem:
    """Minimal stand-in for quiver.System for unit-testing Isotopologue."""
    def __init__(self, hessian, is_linear=False):
        self.hessian = np.asarray(hessian, dtype=float)
        self.number_of_atoms = self.hessian.shape[0] // 3
        self.is_linear = is_linear


# --- mass weighting ----------------------------------------------------------

def test_mass_weighted_hessian():
    # two atoms -> 6x6 hessian; masses 2 and 8 amu
    n = 6
    hessian = np.arange(1, n * n + 1, dtype=float).reshape(n, n)
    hessian = (hessian + hessian.T) / 2.0
    masses = [2.0, 8.0]
    iso = quiver.Isotopologue("t", FakeSystem(hessian), masses)

    m3 = np.repeat(masses, 3)
    expected = hessian / np.sqrt(np.outer(m3, m3))
    assert np.allclose(iso.mw_hessian, expected)


def test_mass_weighting_preserves_symmetry():
    n = 3
    h = np.arange(1, n * n + 1, dtype=float).reshape(n, n)
    h = h + h.T
    iso = quiver.Isotopologue("t", FakeSystem(h), [5.0])
    assert np.allclose(iso.mw_hessian, iso.mw_hessian.T)


# --- frequency calculation + per-object caching ------------------------------

def _make_iso():
    band = np.array([[6, 1, 0, 0, 0, 0],
                     [1, 6, 1, 0, 0, 0],
                     [0, 1, 6, 1, 0, 0],
                     [0, 0, 1, 6, 1, 0],
                     [0, 0, 0, 1, 6, 1],
                     [0, 0, 0, 0, 1, 6]], dtype=float)
    return quiver.Isotopologue("t", FakeSystem(band), [12.0, 12.0])


def test_frequencies_returns_four_tuple():
    small, imag, regular, n_small = _make_iso().calculate_frequencies(
        imag_threshold=50)
    assert isinstance(regular, np.ndarray)
    assert n_small == len(small)


def test_frequencies_second_call_is_cached():
    iso = _make_iso()
    first = iso.calculate_frequencies(imag_threshold=50)
    second = iso.calculate_frequencies(imag_threshold=50)
    # the short-circuit must hand back the very same cached object
    assert first is second


# --- reference isotopologue is diagonalized once, not per substitution --------

def test_reference_diagonalized_once(tutorial):
    cfg = tutorial("gaussian", "claisen_demo.config")
    gs = tutorial("gaussian", "claisen_gs.out")
    ts = tutorial("gaussian", "claisen_ts.out")

    real_eigvalsh = np.linalg.eigvalsh
    calls = []

    def counting_eigvalsh(*a, **k):
        calls.append(1)
        return real_eigvalsh(*a, **k)

    with mock.patch("numpy.linalg.eigvalsh", side_effect=counting_eigvalsh):
        calc = KIE_Calculation(cfg, gs, ts, style="g09")

    n_iso = len(calc.config.isotopologues)
    # cached: 2 reference diagonalizations (gs + ts) + 2 per substitution;
    # the naive recompute-the-reference path would be 4 per substitution
    assert len(calls) == 2 + 2 * n_iso
    assert len(calls) < 4 * n_iso
