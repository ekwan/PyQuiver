"""Tests for the physical-constant and atomic-weight tables in ``constants.py``
(loaded from ``weights.dat``)."""

import pytest

from pyquiver.constants import (DEFAULT_MASSES, REPLACEMENTS, REPLACEMENTS_Z,
                       PHYSICAL_CONSTANTS, LINEARITY_THRESHOLD, DROP_NUM_LINEAR,
                       Element, replacement_mass)


@pytest.mark.parametrize("z,mass", [(1, 1.0078), (6, 12.0), (8, 15.9949)])
def test_default_masses(z, mass):
    assert DEFAULT_MASSES[z] == pytest.approx(mass, abs=1e-3)


@pytest.mark.parametrize("label", ["13C", "2D", "17O"])
def test_known_replacements_present(label):
    assert label in REPLACEMENTS


@pytest.mark.parametrize("label,mass", [("13C", 13.0033), ("2D", 2.0141)])
def test_replacement_masses(label, mass):
    assert REPLACEMENTS[label] == pytest.approx(mass, abs=1e-3)


def test_replacement_mass_resolver():
    assert replacement_mass("13C") == pytest.approx(13.00335, abs=1e-4)  # label
    assert replacement_mass(5000.0) == 5000.0                            # numeric
    for bad in (True, -1.0, "nope"):
        with pytest.raises(ValueError):
            replacement_mass(bad)


def test_replacement_heavier_than_default():
    assert REPLACEMENTS["13C"] > DEFAULT_MASSES[6]


@pytest.mark.parametrize("label,z", [("13C", 6), ("2D", 1), ("17O", 8)])
def test_replacement_atomic_numbers(label, z):
    assert REPLACEMENTS_Z[label] == z


@pytest.mark.parametrize("key", ["h", "c", "kB", "Eh", "a0", "amu", "atb"])
def test_physical_constant_present(key):
    assert key in PHYSICAL_CONSTANTS


@pytest.mark.parametrize("key", ["h", "c", "kB"])
def test_physical_constants_positive(key):
    assert PHYSICAL_CONSTANTS[key] > 0.0


def test_linearity_parameters():
    assert LINEARITY_THRESHOLD > 0.0
    # five modes (3 translation + 2 rotation) are dropped for linear cases
    assert DROP_NUM_LINEAR == 5


# --- Element class -----------------------------------------------------------

def test_element_construction():
    e = Element("carbon", 6, "C", 12.0)
    assert e.full_name == "carbon"
    assert e.symbol == "C"
    assert e.atomic_number == 6
    assert e.default_mass == pytest.approx(12.0)
    assert e.replacements == []


def test_element_add_replacement():
    e = Element("carbon", 6, "C", 12.0)
    e.add_replacement("13C", 13.00335)
    assert ("13C", pytest.approx(13.00335)) in e.replacements


def test_element_str_runs():
    e = Element("carbon", 6, "C", 12.0)
    assert "Carbon" in str(e)
    e.add_replacement("13C", 13.00335)
    assert "13C" in str(e)


def test_element_duplicate_replacement_raises():
    e = Element("carbon", 6, "C", 12.0)
    e.add_replacement("13C", 13.00335)
    with pytest.raises(ValueError):
        e.add_replacement("13C", 13.00335)


@pytest.mark.parametrize("bad", [
    dict(full_name="Carbon", atomic_number=6, symbol="C", default_mass=12.0),   # non-lowercase name
    dict(full_name="carbon", atomic_number=6, symbol="1", default_mass=12.0),   # non-letter symbol
    dict(full_name="carbon", atomic_number=0, symbol="C", default_mass=12.0),   # bad atomic number
    dict(full_name="carbon", atomic_number=6, symbol="CCC", default_mass=12.0), # symbol too long
    dict(full_name="carbon", atomic_number=6, symbol="C", default_mass=999.0),  # mass out of range
])
def test_element_invalid_values_raise(bad):
    with pytest.raises(ValueError):
        Element(**bad)


@pytest.mark.parametrize("symbol", ["&bad", "ABCDE"])
def test_element_add_replacement_validation(symbol):
    e = Element("carbon", 6, "C", 12.0)
    with pytest.raises(ValueError):
        e.add_replacement(symbol, 13.0)


def test_element_add_replacement_bad_mass():
    e = Element("carbon", 6, "C", 12.0)
    with pytest.raises(ValueError):
        e.add_replacement("13C", 999.0)
