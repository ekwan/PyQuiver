"""The optional n_jobs parallelism must not change results.

eigvalsh releases the GIL, so threads parallelize the per-isotopologue work;
this checks that serial and threaded runs agree exactly, and that the
reference isotopologue is still diagonalized only once (pre-warmed before the
fan-out).
"""

from unittest import mock

import numpy as np
import pytest

from pyquiver.kie import KIE_Calculation

CFG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


@pytest.mark.parametrize("n_jobs", [2, 4, -1])
def test_parallel_matches_serial(tutorial, n_jobs):
    serial = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                             style="gaussian", n_jobs=1)
    parallel = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                               style="gaussian", n_jobs=n_jobs)
    assert list(serial.KIES) == list(parallel.KIES)  # order preserved
    for name in serial.KIES:
        assert np.allclose(serial.KIES[name].value, parallel.KIES[name].value)


def test_reference_still_diagonalized_once_in_parallel(tutorial):
    real_eigvalsh = np.linalg.eigvalsh
    calls = []

    def counting(*a, **k):
        calls.append(1)
        return real_eigvalsh(*a, **k)

    with mock.patch("numpy.linalg.eigvalsh", side_effect=counting):
        calc = KIE_Calculation(tutorial(*CFG), tutorial(*GS), tutorial(*TS),
                               style="gaussian", n_jobs=4)

    n_iso = len(calc.config.isotopologues)
    # 2 reference diagonalizations (gs + ts, pre-warmed) + 2 per substitution
    assert len(calls) == 2 + 2 * n_iso
