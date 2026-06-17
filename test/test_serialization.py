"""Round-trip test for the native .qin serialization.

Parsing a g09 file, dumping it to the native ``.qin`` format, and re-parsing
must reproduce the same geometry and Hessian. This exercises
``dump_pyquiver_input_file`` and the ``_make_serial_*`` helpers, and guards
the trailing-comma handling in the lower-triangle Hessian parser.
"""

import shutil

import numpy as np

from pyquiver import quiver


def test_qin_round_trip(tmp_path, tutorial):
    # parse from a copy inside tmp_path so the dumped .qin lands there too
    src = tmp_path / "claisen_gs.out"
    shutil.copy(tutorial("gaussian", "claisen_gs.out"), str(src))
    original = quiver.System(str(src), style="g09")

    original.dump_pyquiver_input_file()  # writes claisen_gs.qin next to src
    qin = tmp_path / "claisen_gs.qin"
    assert qin.exists()

    reparsed = quiver.System(str(qin), style="pyquiver")
    assert reparsed.number_of_atoms == original.number_of_atoms
    assert reparsed.atomic_numbers == original.atomic_numbers
    assert np.allclose(reparsed.hessian, original.hessian, atol=1e-6)
