"""PyQuiver: kinetic and equilibrium isotope effects from Hessian matrices.

PyQuiver computes KIEs/EIEs with the Bigeleisen-Mayer equation from Cartesian
Hessians produced by electronic-structure programs (Gaussian, ORCA, or the
native ``.qin`` format).

Typical programmatic use::

    from pyquiver import KIE_Calculation
    calc = KIE_Calculation("demo.config", "gs.out", "ts.out", style="g09")
    print(calc)
    calc.KIES["C1"].value      # (raw, Wigner, infinite-parabola) KIE
"""

import logging as _logging

from .quiver import System, Isotopologue
from .config import Config
from .kie import KIE_Calculation, KIE
from .results import Results, KIEResult, EIEResult
from .batch import batch, BatchResults
from . import tunneling

# a library should not configure logging itself; attach a no-op handler so
# importing pyquiver never prints "No handlers could be found" and callers
# stay in full control of log output
_logging.getLogger("pyquiver").addHandler(_logging.NullHandler())

try:  # populated from package metadata once installed (PyPI name: pyquiver-kie)
    from importlib.metadata import version as _version
    __version__ = _version("pyquiver-kie")
except Exception:  # pragma: no cover - source checkout that isn't installed
    __version__ = "0+unknown"

__all__ = [
    "System",
    "Isotopologue",
    "Config",
    "KIE_Calculation",
    "KIE",
    "Results",
    "KIEResult",
    "EIEResult",
    "batch",
    "BatchResults",
    "tunneling",
    "__version__",
]
