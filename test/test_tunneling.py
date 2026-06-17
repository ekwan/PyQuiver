"""Tests for the Skodje-Truhlar tunnelling correction (pyquiver.tunneling and
the KIE_Calculation.skodje_truhlar method)."""

import numpy as np
import pytest

from pyquiver.kie import KIE_Calculation
from pyquiver.tunneling import (skodje_truhlar, skodje_truhlar_kappa, wigner,
                                bell, Eh)

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


# --- pure function -----------------------------------------------------------

def test_equal_frequencies_give_unity():
    assert skodje_truhlar(-500.0, -500.0, 300.0, 0.05 * Eh) == pytest.approx(1.0)


def test_rejects_real_frequencies():
    with pytest.raises(ValueError):
        skodje_truhlar(500.0, 500.0, 300.0, 0.05 * Eh)


def test_ratio_matches_kappa_definition():
    barrier = 0.03 * Eh
    expected = (skodje_truhlar_kappa(-400.0, 300.0, barrier)
                / skodje_truhlar_kappa(-500.0, 300.0, barrier))
    assert skodje_truhlar(-500.0, -400.0, 300.0, barrier) == pytest.approx(expected)


def test_kappa_low_temperature_branch():
    # a large imaginary frequency at low temperature gives alpha < beta,
    # exercising the deep-tunnelling form of the equation
    kappa = skodje_truhlar_kappa(-3000.0, 100.0, 0.03 * Eh)
    assert np.isfinite(kappa) and kappa > 0.0
    # the ratio should also evaluate in this regime
    ratio = skodje_truhlar(-3000.0, -2800.0, 100.0, 0.03 * Eh)
    assert np.isfinite(ratio)


def test_module_reexports_wigner_and_bell():
    # all three corrections are reachable from pyquiver.tunneling
    assert wigner(-500.0, -500.0, 300.0) == pytest.approx(1.0)
    assert bell(-500.0, -500.0, 300.0) == pytest.approx(1.0)


# --- high-level method -------------------------------------------------------

@pytest.fixture
def calc(tutorial):
    return KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")


def test_method_returns_all_isotopologues(calc):
    st = calc.skodje_truhlar(reactant_energy=0.0, product_energy=0.0,
                             ts_energy=0.02)
    assert set(st) == set(calc.KIES)
    assert all(np.isfinite(v) for v in st.values())


def test_hartree_and_joule_units_agree(calc):
    a = calc.skodje_truhlar(0.0, 0.0, 0.02, unit="hartree")
    b = calc.skodje_truhlar(0.0, 0.0, 0.02 * Eh, unit="J")
    assert a == pytest.approx(b)


def test_bad_unit_rejected(calc):
    with pytest.raises(ValueError):
        calc.skodje_truhlar(0.0, 0.0, 0.02, unit="kcal")


def test_method_without_reference(tutorial):
    # reference 'none' -> absolute ST-corrected KIEs (exercises the no-ref path)
    from pyquiver.config import Config
    from pyquiver import System
    cfg = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                           temperature=393, scaling=0.9614, imag_threshold=50)
    gs = System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    calc = KIE_Calculation(cfg, gs, ts)
    st = calc.skodje_truhlar(0.0, 0.0, 0.02)
    assert np.isfinite(st["C1"])


def test_eie_calculation_rejected(tutorial):
    eie = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*GS),
                          style="gaussian")
    with pytest.raises(ValueError):
        eie.skodje_truhlar(0.0, 0.0, 0.02)
