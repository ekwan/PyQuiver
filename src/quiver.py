"""Deprecated shim. Import ``pyquiver.quiver`` (or use the ``pyquiver``
command / ``python -m pyquiver``) instead of running ``src/quiver.py``.

Kept so that pre-packaging usage keeps working:
    python src/quiver.py CONFIG GS TS
    from quiver import System, Isotopologue
"""
import os as _os, sys as _sys, warnings as _warnings

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

if __name__ == "__main__":
    _warnings.warn("Use the 'pyquiver' command or 'python -m pyquiver' "
                   "instead of 'python src/quiver.py'.",
                   DeprecationWarning, stacklevel=2)
    from pyquiver.cli import main
    main()
else:
    _warnings.warn("'quiver' is deprecated; import 'pyquiver.quiver' instead.",
                   DeprecationWarning, stacklevel=2)
    import pyquiver.quiver as _mod
    _sys.modules[__name__] = _mod
