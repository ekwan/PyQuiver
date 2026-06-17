"""Tests for the KIE/EIE orchestration in ``kie.py`` (the KIE class and the
KIE_Calculation driver), beyond the pure physics functions."""

import numpy as np
import pytest

from pyquiver.kie import KIE_Calculation
from pyquiver.config import Config
from pyquiver import quiver

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


@pytest.fixture
def kie(tutorial):
    """Run a KIE_Calculation; pass which tutorial files to use as gs and ts."""
    def _run(gs=GS, ts=TS):
        return KIE_Calculation(tutorial(*CFG), tutorial(*gs), tutorial(*ts),
                               style="g09")
    return _run


def test_transition_state_is_kie(kie):
    # a real transition state has an imaginary mode -> KIE (eie_flag 0), and
    # each value is a 3-vector (raw, Wigner, inverted parabola)
    calc = kie()
    assert calc.eie_flag == 0
    assert len(calc.KIES["C1"].value) == 3


def test_no_imaginary_mode_is_eie(kie):
    # using the ground state as both endpoints -> no imaginary mode -> EIE
    calc = kie(gs=GS, ts=GS)
    assert calc.eie_flag == 1
    assert np.ndim(calc.KIES["C1"].value) == 0


def test_eie_of_identical_endpoints_is_unity(kie):
    calc = kie(gs=GS, ts=GS)
    assert float(calc.KIES["C1"].value) == pytest.approx(1.0, abs=1e-6)


def test_reference_isotopologue_kept_absolute(kie):
    # the config references everything to C5; the reference entry itself is
    # left as its absolute value (~1.0019 for the raw column), not divided
    calc = kie()
    assert calc.config.reference_isotopologue == "C5"
    assert float(calc.KIES["C5"].value[0]) == pytest.approx(1.001895, abs=1e-5)


@pytest.mark.parametrize("endpoints", [(GS, TS), (GS, GS)])
def test_str_runs(kie, endpoints):
    assert "Isotopologue" in str(kie(*endpoints))


def test_get_row_shapes(kie):
    title_row, row, eie_p = kie().get_row()
    assert eie_p == 0
    # the reference (C5) is excluded from the output columns
    assert "C5" not in title_row
    # one comma-separated value per reported isotopologue
    n_titles = len([t for t in title_row.split(",") if t])
    n_values = len([v for v in row.split(",") if v])
    assert n_titles == n_values


def test_accepts_prebuilt_objects(tutorial):
    # KIE_Calculation accepts Config and System objects, not just file paths;
    # the result must match the path-based calculation
    cfg = Config(tutorial(*CFG))
    gs = quiver.System(tutorial(*GS), style="g09")
    ts = quiver.System(tutorial(*TS), style="g09")
    from_objects = KIE_Calculation(cfg, gs, ts)
    from_paths = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                                 style="g09")
    assert [float(v) for v in from_objects.KIES["C1"].value] == \
        pytest.approx([float(v) for v in from_paths.KIES["C1"].value])


def test_rejects_bad_argument_types(tutorial):
    with pytest.raises(TypeError):
        KIE_Calculation(123, tutorial(*GS), tutorial(*TS), style="g09")


def test_apply_reference_divides(kie):
    calc = kie()
    # apply_reference divides this KIE's value in place by the reference's;
    # dividing by itself yields ones in every column
    result = calc.KIES["C1"].apply_reference(calc.KIES["C1"])
    assert np.allclose(result, 1.0)
