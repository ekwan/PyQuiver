"""Tests for the structured result objects (pyquiver.results)."""

import pytest

from pyquiver.kie import KIE_Calculation
from pyquiver.results import KIEResult, EIEResult

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


@pytest.fixture
def calc(tutorial):
    return KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                           style="gaussian")


def test_named_access(calc):
    r = calc.results["C1"]
    assert isinstance(r, KIEResult)
    assert r.name == "C1"
    # the three columns match the positional .value vector
    assert (r.uncorrected, r.wigner, r.inverted_parabola) == \
        pytest.approx(tuple(float(v) for v in calc.KIES["C1"].value))


def test_to_dict_kie(calc):
    d = calc.to_dict()
    assert set(d) == set(calc.KIES)
    assert set(d["C1"]) == {"uncorrected", "wigner", "inverted_parabola"}
    assert d["C1"]["uncorrected"] == pytest.approx(1.010830, abs=1e-5)


def test_to_csv_roundtrip(calc, tmp_path):
    path = tmp_path / "out.csv"
    text = calc.to_csv(str(path))
    assert path.read_text() == text
    lines = [l for l in text.splitlines() if l.strip()]
    assert lines[0] == "name,uncorrected,wigner,inverted_parabola"
    # one data row per isotopologue
    assert len(lines) == 1 + len(calc.KIES)


def test_to_dataframe(calc):
    pandas = pytest.importorskip("pandas")
    df = calc.to_dataframe()
    assert list(df.columns) == ["name", "uncorrected", "wigner", "inverted_parabola"]
    assert len(df) == len(calc.KIES)
    row = df[df["name"] == "C1"].iloc[0]
    assert row["wigner"] == pytest.approx(1.012486, abs=1e-5)


def test_eie_results_are_scalar(tutorial):
    # gs as both endpoints -> EIE -> single value column
    calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*GS),
                           style="gaussian")
    assert calc.results.is_eie
    r = calc.results["C1"]
    assert isinstance(r, EIEResult)
    assert calc.to_dict()["C1"] == pytest.approx(1.0, abs=1e-6)
    assert calc.results.columns == ["name", "value"]
