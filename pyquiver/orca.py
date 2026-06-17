"""Backward-compatibility re-export.

The ORCA parser now lives in :mod:`pyquiver.parsers.orca`. This module keeps
``pyquiver.orca.parse_orca_output`` importable for older code.
"""

from .parsers.orca import parse_orca_output, parse

__all__ = ["parse_orca_output", "parse"]
