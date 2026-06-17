"""Error-path and less-common-branch coverage."""

import logging

import numpy as np
import pytest

from pyquiver import quiver
from pyquiver.config import Config
from pyquiver.kie import KIE_Calculation
from pyquiver.kie import calculate_rpfr
from pyquiver.parsers.orca import parse_orca_output

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


# --- Results container -------------------------------------------------------

def test_results_iter_len_and_missing(tutorial):
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")
    res = calc.results
    assert len(res) == len(list(res))
    # the reference isotopologue is excluded from the structured results
    assert {r.name for r in res} == set(calc.KIES) - {calc.config.reference_isotopologue}
    with pytest.raises(KeyError):
        res["does-not-exist"]


def test_to_dataframe_without_pandas(monkeypatch, tutorial):
    # simulate pandas being unavailable
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "pandas":
            raise ImportError("no pandas")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")
    with pytest.raises(ImportError):
        calc.to_dataframe()


# --- KIE_Calculation argument validation -------------------------------------

def test_bad_config_type(tutorial):
    with pytest.raises(TypeError):
        KIE_Calculation(123, tutorial(*GS), tutorial(*TS))


def test_bad_gs_type(tutorial):
    with pytest.raises(TypeError):
        KIE_Calculation(tutorial(*CFG), 123, tutorial(*TS), style="gaussian")


def test_bad_ts_type(tutorial):
    with pytest.raises(TypeError):
        KIE_Calculation(tutorial(*CFG), tutorial(*GS), 123, style="gaussian")


def test_frequency_threshold_warning(tutorial, caplog):
    # a config carrying the deprecated frequency_threshold warns when used
    cfg = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                           temperature=393, scaling=0.9614, imag_threshold=50)
    cfg.frequency_threshold = 50.0
    gs = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = quiver.System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    with caplog.at_level(logging.WARNING, logger="pyquiver"):
        KIE_Calculation(cfg, gs, ts)
    assert any("frequency_threshold" in m for m in caplog.messages)


# --- get_row with tunnelling + debug logging branches ------------------------

def test_get_row_report_tunnelling(tutorial):
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")
    title, row, eie = calc.get_row(report_tunnelling=True)
    assert "_uncorr" in title and "_inf_para" in title
    assert eie == 0


def test_debug_logging_paths(tutorial, caplog):
    # exercise the logger.isEnabledFor(DEBUG) branches in calculate_rpfr /
    # partition_components and the per-mode dump
    with caplog.at_level(logging.DEBUG, logger="pyquiver"):
        KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                        style="gaussian")
    assert any("MODE" in m or "Factor" in m for m in caplog.messages)


# --- mass override isotopologue path -----------------------------------------

def test_mass_override_isotopologue(tutorial):
    cfg = Config.from_dict(
        isotopologues={"C1": [(1, 1, "13C")], "C4": [(4, 4, "13C")]},
        temperature=393, scaling=0.9614, imag_threshold=50,
        reference_isotopologue="C1", mass_override_isotopologue="C4")
    gs = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = quiver.System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    calc = KIE_Calculation(cfg, gs, ts)
    assert "C1" in calc.KIES
    # get_row and __str__ both drop the reference and mass-override columns
    title, row, eie = calc.get_row()
    assert "C4" not in title and "C1" not in title
    assert "Isotopologue" in str(calc)


def test_str_with_and_without_reference(tutorial):
    # with a reference (C5) -> "referenced to" section
    ref = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                          style="gaussian")
    assert "referenced to isotopologue C5" in str(ref)

    # reference "none" -> "Absolute KIEs" section
    cfg = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                           temperature=393, scaling=0.9614, imag_threshold=50)
    gs = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = quiver.System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    assert "Absolute KIEs" in str(KIE_Calculation(cfg, gs, ts))


def test_get_row_eie(tutorial):
    # EIE (gs as both endpoints) -> single value column in get_row
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*GS),
                           style="gaussian")
    title, row, eie = calc.get_row()
    assert eie == 1
    assert row  # EIE values rendered


def test_build_default_masses_missing_mass(tutorial):
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")

    class Stub:
        number_of_atoms = 1
        atomic_numbers = [200]   # not present in weights.dat
        filename = "stub"

    with pytest.raises(ValueError):
        calc.build_default_masses(Stub())


# --- Isotopologue / frequency branches ---------------------------------------

class _FakeSystem:
    def __init__(self, hessian, is_linear=False):
        self.hessian = np.asarray(hessian, dtype=float)
        self.number_of_atoms = self.hessian.shape[0] // 3
        self.is_linear = is_linear


def test_isotopologue_str():
    iso = quiver.Isotopologue("x", _FakeSystem(np.eye(3)), [12.0])
    assert "Isotopologue: x" in str(iso)


def test_calculate_frequencies_unknown_method():
    iso = quiver.Isotopologue("x", _FakeSystem(np.eye(6)), [12.0, 12.0])
    with pytest.raises(ValueError):
        iso.calculate_frequencies(imag_threshold=50, method="bogus")


def test_calculate_frequencies_linear_branch():
    iso = quiver.Isotopologue("x", _FakeSystem(np.eye(6), is_linear=True),
                              [12.0, 12.0])
    result = iso.calculate_frequencies(imag_threshold=50)
    assert len(result) == 4


def test_multiple_imaginary_warning(caplog):
    # a Hessian with two strongly negative eigenvalues yields two imaginary
    # frequencies, which logs a warning
    hess = np.diag([-1.0, -1.0, 1.0, 1.0, 1.0, 1.0])
    iso = quiver.Isotopologue("x", _FakeSystem(hess), [12.0, 12.0])
    with caplog.at_level(logging.WARNING, logger="pyquiver"):
        iso.calculate_frequencies(imag_threshold=50)
    assert any("imaginary" in m for m in caplog.messages)


# --- native (.qin) parser errors ---------------------------------------------

def test_native_bad_atom_count(tmp_path):
    p = tmp_path / "bad.qin"
    p.write_text("notanumber\n")
    with pytest.raises(ValueError):
        quiver.System(str(p), style="native")


def test_native_bad_geometry_line(tmp_path):
    p = tmp_path / "bad.qin"
    p.write_text("1\n0,6,0.0\n0.0,\n")   # geometry line has too few fields
    with pytest.raises(ValueError):
        quiver.System(str(p), style="native")


# --- ORCA parser errors ------------------------------------------------------

def test_orca_missing_atoms():
    with pytest.raises(ValueError):
        parse_orca_output("nothing useful here")


def test_orca_missing_hessian():
    data = "$atoms\n1\n C 1.0 0.0 0.0 0.0\n"
    with pytest.raises(ValueError):
        parse_orca_output(data)


# --- calculate_rpfr imaginary-frequency-count mismatch -----------------------

class _FreqStub:
    def __init__(self, small, imag, freqs):
        self._r = (small, imag, np.array(freqs, dtype=float), len(small))

    def calculate_frequencies(self, imag_threshold, scaling=1.0):
        return self._r


def test_calculate_rpfr_return_order(tutorial):
    # contract: (rpfr, imag_ratios, heavy_freqs, light_freqs) -- heavy at index
    # 2, light at index 3. The heavier isotopologue has the lower frequencies.
    from pyquiver import System, Isotopologue
    from pyquiver.kie import calculate_rpfr
    from pyquiver.constants import DEFAULT_MASSES, REPLACEMENTS
    gs = System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    light_masses = [DEFAULT_MASSES[z] for z in gs.atomic_numbers]
    heavy_masses = list(light_masses)
    heavy_masses[0] = REPLACEMENTS["13C"]
    light = Isotopologue("light", gs, light_masses)
    heavy = Isotopologue("heavy", gs, heavy_masses)
    rpfr, imag, heavy_freqs, light_freqs = calculate_rpfr((light, heavy), 50.0, 1.0, 393)
    assert (heavy_freqs <= light_freqs + 1e-9).all()   # heavy never higher


def test_calculate_rpfr_imaginary_mismatch_warns(caplog):
    # light has one imaginary mode, heavy has none -> warn and drop it
    light = _FreqStub([], [-500.0], [1000.0, 1100.0])
    heavy = _FreqStub([], [], [1000.0, 1100.0])
    with caplog.at_level(logging.WARNING, logger="pyquiver"):
        rpfr, imag_ratios, _, _ = calculate_rpfr((light, heavy), 50, 1.0, 300.0)
    assert any("imaginary" in m for m in caplog.messages)
    assert imag_ratios is None   # imaginary mode was dropped
