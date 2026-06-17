"""Deprecated shim. Import `pyquiver.config` instead of `config`.

Kept so that pre-packaging scripts (`from config import ...`) keep working.
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_warnings.warn(
    "'config' is deprecated; import 'pyquiver.config' instead.",
    DeprecationWarning, stacklevel=2)

import pyquiver.config as _mod
_sys.modules[__name__] = _mod
