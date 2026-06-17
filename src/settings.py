"""Deprecated shim. Import `pyquiver.settings` instead of `settings`.

Kept so that pre-packaging scripts (`from settings import ...`) keep working.
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_warnings.warn(
    "'settings' is deprecated; import 'pyquiver.settings' instead.",
    DeprecationWarning, stacklevel=2)

import pyquiver.settings as _mod
_sys.modules[__name__] = _mod
