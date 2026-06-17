"""Deprecated shim. Import `pyquiver.orca` instead of `orca`.

Kept so that pre-packaging scripts (`from orca import ...`) keep working.
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_warnings.warn(
    "'orca' is deprecated; import 'pyquiver.orca' instead.",
    DeprecationWarning, stacklevel=2)

import pyquiver.orca as _mod
_sys.modules[__name__] = _mod
