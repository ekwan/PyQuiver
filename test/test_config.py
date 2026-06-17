"""Tests for the configuration-file parser (``config.Config``)."""

import pytest

from pyquiver.config import Config
from pyquiver import quiver

MINIMAL = """\
scaling 0.96
temperature 300
imag_threshold 50
mass_override_isotopologue default
reference_isotopologue none
isotopologue C1 1 1 13C
"""


@pytest.fixture
def make_config(tmp_path):
    """Return a helper that writes config text to a temp file and parses it."""
    counter = {"n": 0}

    def _make(text):
        counter["n"] += 1
        path = tmp_path / ("c%d.config" % counter["n"])
        path.write_text(text)
        return Config(str(path))
    return _make


def test_scalar_fields(make_config):
    c = make_config(MINIMAL)
    assert c.temperature == 300.0
    assert c.scaling == 0.96
    assert c.imag_threshold == 50.0
    assert c.reference_isotopologue == "none"


def test_single_isotopologue(make_config):
    c = make_config(MINIMAL)
    assert list(c.isotopologues) == ["C1"]
    assert c.isotopologues["C1"] == [(1, 1, "13C")]


def test_deprecated_isotopomer_alias(make_config):
    text = MINIMAL.replace("isotopologue C1", "isotopomer C1") \
                  .replace("reference_isotopologue none",
                           "reference_isotopomer none")
    assert "C1" in make_config(text).isotopologues


def test_multiple_replacements_one_isotopologue(make_config):
    c = make_config(MINIMAL + "isotopologue C1 2 2 13C\n")
    assert c.isotopologues["C1"] == [(1, 1, "13C"), (2, 2, "13C")]


def test_comments_and_blank_lines_ignored(make_config):
    c = make_config("# a comment\n\n" + MINIMAL + "  # trailing\n")
    assert c.scaling == 0.96


@pytest.mark.parametrize("mutation", [
    "isotopologue default 2 2 13C\n",          # reserved name
    "isotopologue C2 2 2 99Z\n",               # invalid replacement
])
def test_bad_isotopologue_rejected(make_config, mutation):
    with pytest.raises(ValueError):
        make_config(MINIMAL + mutation)


@pytest.mark.parametrize("old,new", [
    ("scaling 0.96", "scaling 2.0"),           # out of range
    ("temperature 300", "temperature -5"),     # negative
    ("imag_threshold 50", "imag_threshold 500"),  # too high
])
def test_out_of_range_fields_rejected(make_config, old, new):
    with pytest.raises(ValueError):
        make_config(MINIMAL.replace(old, new))


@pytest.mark.parametrize("prefix", ["temperature", "isotopologue"])
def test_missing_required_lines_rejected(make_config, prefix):
    text = "\n".join(l for l in MINIMAL.splitlines()
                     if not l.startswith(prefix)) + "\n"
    with pytest.raises(ValueError):
        make_config(text)


def test_str_runs(make_config):
    assert "Config file" in str(make_config(MINIMAL))


def test_mass_override_isotopomer_alias(make_config):
    text = MINIMAL + "mass_override_isotopomer C1\n"
    assert make_config(text).mass_override_isotopologue == "C1"


@pytest.mark.parametrize("bad_line", [
    "isotopologue C1 1 1\n",        # too few fields for an isotopologue
    "isotopologue C2 0 1 13C\n",    # atom number < 1
    "bogus_field value\n",          # unknown 2-field config key
    "one two three\n",              # wrong number of fields entirely
])
def test_malformed_lines_rejected(make_config, bad_line):
    with pytest.raises(ValueError):
        make_config(MINIMAL + bad_line)


def test_frequency_threshold_warned_and_shown(make_config, caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="pyquiver"):
        c = make_config(MINIMAL + "frequency_threshold 50\n")
    assert any("frequency_threshold" in m for m in caplog.messages)
    assert "Frequency threshold" in str(c)


def test_str_reference_not_found_raises(make_config):
    c = make_config(MINIMAL.replace("reference_isotopologue none",
                                    "reference_isotopologue ZZ"))
    with pytest.raises(ValueError):
        str(c)


def test_from_dict_numeric_mass(tutorial):
    # a bare number is used directly as an unusual mass (no weights.dat label)
    from pyquiver.kie import KIE_Calculation
    from pyquiver import System
    cfg = Config.from_dict(isotopologues={"heavyC4": [(4, 4, 5000.0)]},
                           temperature=393, scaling=0.9614, imag_threshold=50)
    assert cfg.isotopologues["heavyC4"] == [(4, 4, 5000.0)]   # numeric type preserved
    gs = System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    calc = KIE_Calculation(cfg, gs, ts)
    assert calc.results["heavyC4"].uncorrected > 1.5   # very heavy -> large KIE


def test_from_dict_negative_mass_rejected():
    with pytest.raises(ValueError):
        Config.from_dict(isotopologues={"x": [(1, 1, -3.0)]},
                         temperature=300, scaling=0.96, imag_threshold=50)


def test_config_file_numeric_mass(make_config):
    c = make_config(MINIMAL + "isotopologue heavyC4 4 4 5000.0\n")
    assert c.isotopologues["heavyC4"] == [(4, 4, 5000.0)]


def test_check_skips_element_match_for_numeric_mass(make_config, claisen_systems):
    # a custom numeric mass carries no element identity, so check() must not
    # try to validate the element and must not raise
    gs, ts = claisen_systems
    c = make_config(MINIMAL + "isotopologue heavyC4 4 4 5000.0\n")
    c.check(gs, ts)   # must not raise


def test_from_dict_atom_index_validated():
    with pytest.raises(ValueError):
        Config.from_dict(isotopologues={"C1": [(0, 1, "13C")]},
                         temperature=300, scaling=0.96, imag_threshold=50)


def test_check_verbose_logs(make_config, claisen_systems, caplog):
    import logging
    gs, ts = claisen_systems
    with caplog.at_level(logging.INFO, logger="pyquiver"):
        make_config(MINIMAL).check(gs, ts, verbose=True)
    assert any("makes sense" in m for m in caplog.messages)


def test_check_ts_mismatch_rejected(make_config):
    # gs atom matches the replacement element but the ts atom does not
    class Stub:
        def __init__(self, znums):
            self.atomic_numbers = znums
            self.filename = "stub"
    c = make_config(MINIMAL)  # isotopologue C1: gs atom 1 -> 13C (Z=6)
    gs = Stub([6])   # carbon, matches
    ts = Stub([7])   # nitrogen, mismatch
    with pytest.raises(ValueError):
        c.check(gs, ts)


def test_check_valid_config_passes(make_config, claisen_systems):
    gs, ts = claisen_systems
    make_config(MINIMAL).check(gs, ts)  # must not raise


def test_from_dict_builds_without_file():
    c = Config.from_dict(
        isotopologues={"C1": [(1, 1, "13C")]},
        temperature=300, scaling=0.96, imag_threshold=50)
    assert c.temperature == 300.0
    assert c.scaling == 0.96
    assert c.reference_isotopologue == "none"
    assert c.isotopologues["C1"] == [(1, 1, "13C")]


def test_from_dict_matches_file(make_config):
    from_file = make_config(MINIMAL)
    from_dict = Config.from_dict(
        isotopologues={"C1": [(1, 1, "13C")]},
        temperature=300, scaling=0.96, imag_threshold=50)
    assert from_dict.isotopologues == from_file.isotopologues
    assert from_dict.temperature == from_file.temperature
    assert from_dict.scaling == from_file.scaling


@pytest.mark.parametrize("kwargs", [
    dict(isotopologues={"C1": [(1, 1, "99Z")]}, temperature=300, scaling=0.96, imag_threshold=50),
    dict(isotopologues={"default": [(1, 1, "13C")]}, temperature=300, scaling=0.96, imag_threshold=50),
    dict(isotopologues={"C1": [(1, 1, "13C")]}, temperature=300, scaling=2.0, imag_threshold=50),
    dict(isotopologues={}, temperature=300, scaling=0.96, imag_threshold=50),
])
def test_from_dict_validates(kwargs):
    with pytest.raises(ValueError):
        Config.from_dict(**kwargs)


def test_from_dict_usable_in_calculation(tutorial):
    # a programmatic Config should drive a real calculation
    from pyquiver.kie import KIE_Calculation
    cfg = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                           temperature=393, scaling=0.9614, imag_threshold=50)
    gs = quiver.System(tutorial("gaussian", "claisen_gs.out"), style="gaussian")
    ts = quiver.System(tutorial("gaussian", "claisen_ts.out"), style="gaussian")
    calc = KIE_Calculation(cfg, gs, ts)
    # reference defaults to "none", so this is the absolute (un-referenced) KIE
    assert calc.results["C1"].uncorrected == pytest.approx(1.0127464, abs=1e-5)


def test_check_element_mismatch_rejected(make_config, claisen_systems):
    # atom 1 of the Claisen system is carbon; labelling it 2D (hydrogen) is
    # an element mismatch that check() must reject
    gs, ts = claisen_systems
    c = make_config(MINIMAL.replace("isotopologue C1 1 1 13C",
                                    "isotopologue X 1 1 2D"))
    with pytest.raises(ValueError):
        c.check(gs, ts)
