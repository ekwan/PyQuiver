"""Tests for the command-line entry point (pyquiver.cli)."""

import pytest

from pyquiver import cli

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


def test_main_runs_and_prints(tutorial, capsys):
    calc = cli.main([tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                     "-s", "gaussian"])
    out = capsys.readouterr().out
    assert "Isotopologue" in out
    assert "C1" in calc.KIES


@pytest.mark.parametrize("extra", [["-v"], ["-vv"], ["-j", "2"]])
def test_main_flag_variants(tutorial, extra):
    calc = cli.main([tutorial(*CFG), tutorial(*GS), tutorial(*TS)] + extra)
    assert calc.eie_flag == 0


def test_main_bad_style_raises(tutorial):
    with pytest.raises(ValueError):
        cli.main([tutorial(*CFG), tutorial(*GS), tutorial(*TS), "-s", "nope"])


def test_main_is_deprecated(tutorial, capsys):
    with pytest.warns(DeprecationWarning):
        cli.main([tutorial(*CFG), tutorial(*GS), tutorial(*TS)])
    assert "deprecated" in capsys.readouterr().err
