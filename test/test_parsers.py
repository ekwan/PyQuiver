"""Tests for the input-file parsers (``quiver.System``).

The Claisen tutorial system has 14 atoms, so every parser should produce a
42x42 symmetric Hessian.
"""

import numpy as np
import pytest

from pyquiver import quiver
from pyquiver import parsers

N_ATOMS = 14


@pytest.mark.parametrize("style", ["gaussian", "g16", "g09"])
def test_gaussian_style_aliases(style, tutorial):
    # gaussian is the preferred name; g16/g09 are accepted aliases
    system = quiver.System(tutorial("gaussian", "claisen_gs.out"), style=style)
    assert system.number_of_atoms == N_ATOMS


@pytest.mark.parametrize("style", ["pyquiver", "native", "qin"])
def test_native_style_aliases(style, tutorial):
    system = quiver.System(tutorial("pyquiver", "claisen_gs.qin"), style=style)
    assert system.number_of_atoms == N_ATOMS


def test_supported_styles_listed():
    styles = parsers.supported_styles()
    for expected in ("gaussian", "orca", "native"):
        assert expected in styles


@pytest.mark.parametrize("style,filename", [
    ("g09", "gaussian/claisen_gs.out"),
    ("g09", "gaussian/claisen_ts.out"),
    ("orca", "orca/claisen_gs_freq.hess"),
    ("orca", "orca/claisen_ts_freq.hess"),
    ("pyquiver", "pyquiver/claisen_gs.qin"),
    ("pyquiver", "pyquiver/claisen_ts.qin"),
])
def test_parsers_produce_symmetric_hessian(style, filename, tutorial):
    system = quiver.System(tutorial(*filename.split("/")), style=style)
    assert system.number_of_atoms == N_ATOMS
    assert system.hessian.shape == (3 * N_ATOMS, 3 * N_ATOMS)
    assert np.allclose(system.hessian, system.hessian.T), \
        "parsed Hessian is not symmetric"


def test_invalid_style_rejected(tutorial):
    with pytest.raises(ValueError):
        quiver.System(tutorial("gaussian", "claisen_gs.out"),
                      style="not_a_real_style")


def test_g09_and_orca_agree_on_size(tutorial):
    g = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="g09")
    o = quiver.System(tutorial("orca", "claisen_gs_freq.hess"), style="orca")
    assert g.number_of_atoms == o.number_of_atoms


def test_pyquiver_round_trips_against_g09(tutorial):
    # the .qin files were generated from the g09 output, so the parsed
    # Hessians should match what the g09 parser produces
    qin = quiver.System(tutorial("pyquiver", "claisen_gs.qin"),
                        style="pyquiver")
    g09 = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="g09")
    assert np.allclose(qin.hessian, g09.hessian, atol=1e-6)
