"""Deprecated shim. Import `pyquiver.kie` instead of `kie`.

Kept so that pre-packaging scripts (`from kie import ...`) keep working.
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_warnings.warn(
    "'kie' is deprecated; import 'pyquiver.kie' instead.",
    DeprecationWarning, stacklevel=2)

import pyquiver.kie as _mod
_sys.modules[__name__] = _mod
