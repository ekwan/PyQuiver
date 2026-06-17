"""Error and fallback paths in the Gaussian parser, using crafted minimal
output files."""

import pytest

from pyquiver import quiver
from pyquiver.parsers import gaussian


def _write(tmp_path, text, name="job.out"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_tail_handles_tiny_file(tmp_path):
    # a file too short to seek backwards exercises the OSError fallback in _tail
    p = tmp_path / "tiny.out"
    p.write_text("x")
    assert gaussian._tail(str(p)) == "x"


def test_terminated_in_error(tmp_path):
    path = _write(tmp_path, "some output\njob died\n")
    with pytest.raises(ValueError):
        quiver.System(path, style="gaussian")


def test_missing_verbose_flag(tmp_path):
    path = _write(tmp_path, "output line\nNormal termination of Gaussian\n")
    with pytest.raises(ValueError):
        quiver.System(path, style="gaussian")


def test_missing_natoms(tmp_path):
    path = _write(tmp_path, " #p \nNormal termination\n")
    with pytest.raises(ValueError):
        quiver.System(path, style="gaussian")


def test_missing_geometry(tmp_path):
    path = _write(tmp_path, " #p \nNAtoms=    1 \nNormal termination\n")
    with pytest.raises(ValueError):
        quiver.System(path, style="gaussian")


# a standard-orientation geometry but no frequency archive -> ValueError, after
# exercising the geometry-table parsing
STD_ORIENT = (
    " #p \n"
    "NAtoms=    1 \n"
    "Standard orientation\n"
    "    1    6    0    0.000000    0.000000    0.000000\n"
    " Rotational constants (GHZ)\n"
    "Normal termination\n"
)


def test_no_frequency_archive(tmp_path):
    with pytest.raises(ValueError):
        quiver.System(_write(tmp_path, STD_ORIENT), style="gaussian")


# input-orientation fallbacks (no standard orientation present)
INPUT_ORIENT_DISTANCE = (
    " #p \n"
    "NAtoms=    1 \n"
    "Input orientation\n"
    "    1    6    0    0.000000    0.000000    0.000000\n"
    " Distance matrix\n"
    "Normal termination\n"
)

INPUT_ORIENT_ROT = (
    " #p \n"
    "NAtoms=    1 \n"
    "Input orientation\n"
    "    1    6    0    0.000000    0.000000    0.000000\n"
    " Rotational constants (GHZ)\n"
    "Normal termination\n"
)


@pytest.mark.parametrize("text", [INPUT_ORIENT_DISTANCE, INPUT_ORIENT_ROT])
def test_input_orientation_fallbacks(tmp_path, text):
    # geometry is found via input orientation, then the missing archive raises
    with pytest.raises(ValueError):
        quiver.System(_write(tmp_path, text), style="gaussian")
