"""End-to-end regression tests.

These run a full ``KIE_Calculation`` on the tutorial Claisen rearrangement
inputs and lock in the reported KIE values. They are the safety net for any
refactoring (performance work, vectorization, etc.): if the numbers move, a
test breaks.

Golden values were captured from the current implementation and rounded to
six decimals; the printed output is only reported to four decimals, so the
1e-5 tolerance is well within what anyone reads off a result.
"""

import pytest

from pyquiver.kie import KIE_Calculation

# Each entry: style -> (config, gs, ts files) and the expected
# (raw, Wigner, infinite-parabola) KIE per isotopologue. Non-reference
# isotopologues are already divided by the C5 reference.
GOLDEN = {
    "g09": {
        "files": ("gaussian/claisen_demo.config",
                  "gaussian/claisen_gs.out",
                  "gaussian/claisen_ts.out"),
        "values": {
            "C1":  (1.010830, 1.012486, 1.012778),
            "C2":  (1.000021, 1.000043, 1.000047),
            "O3":  (1.016914, 1.018397, 1.018659),
            "C4":  (1.027808, 1.030531, 1.031010),
            "C5":  (1.001895, 1.001933, 1.001940),
            "C6":  (1.012881, 1.014584, 1.014886),
            "H/D": (0.953030, 0.954466, 0.954720),
        },
    },
    "orca": {
        "files": ("orca/claisen_demo.config",
                  "orca/claisen_gs_freq.hess",
                  "orca/claisen_ts_freq.hess"),
        "values": {
            "C1":  (1.011186, 1.012815, 1.013100),
            "C2":  (1.000320, 1.000343, 1.000347),
            "O3":  (1.017148, 1.018604, 1.018860),
            "C4":  (1.028062, 1.030778, 1.031251),
            "C5":  (1.001825, 1.001862, 1.001868),
            "C6":  (1.012909, 1.014593, 1.014888),
            "H/D": (0.947056, 0.948577, 0.948843),
        },
    },
    # The native .qin inputs are generated from the g09 output, so the KIEs
    # are identical to the g09 case (this also exercises the .qin parser).
    "pyquiver": {
        "files": ("pyquiver/claisen_demo.config",
                  "pyquiver/claisen_gs.qin",
                  "pyquiver/claisen_ts.qin"),
        "values": {
            "C1":  (1.010830, 1.012486, 1.012778),
            "C2":  (1.000021, 1.000043, 1.000047),
            "O3":  (1.016914, 1.018397, 1.018659),
            "C4":  (1.027808, 1.030531, 1.031010),
            "C5":  (1.001895, 1.001933, 1.001940),
            "C6":  (1.012881, 1.014584, 1.014886),
            "H/D": (0.953030, 0.954466, 0.954720),
        },
    },
}


@pytest.mark.parametrize("style", list(GOLDEN))
def test_kies_match_golden(style, tutorial):
    cfg, gs, ts = (tutorial(*p.split("/")) for p in GOLDEN[style]["files"])
    calc = KIE_Calculation(cfg, gs, ts, style=style)

    expected = GOLDEN[style]["values"]
    assert set(calc.KIES) == set(expected), "set of isotopologues changed"

    for name, exp in expected.items():
        got = [float(v) for v in calc.KIES[name].value]
        assert got == pytest.approx(exp, abs=1e-5), \
            "%s style, isotopologue %s: %s != %s" % (style, name, got, exp)
