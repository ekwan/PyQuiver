"""Deprecated shim. Import `pyquiver.constants` instead of `constants`.

Kept so that pre-packaging scripts (`from constants import ...`) keep working.
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_warnings.warn(
    "'constants' is deprecated; import 'pyquiver.constants' instead.",
    DeprecationWarning, stacklevel=2)

import pyquiver.constants as _mod
_sys.modules[__name__] = _mod
