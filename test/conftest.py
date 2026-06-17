"""Shared pytest configuration and fixtures for the PyQuiver test suite.

The ``pyquiver`` package lives at the repository root, so adding the root to
``sys.path`` makes ``import pyquiver`` resolve directly from the source tree
without requiring an install. pytest captures stdout by default, so the
package's chatty prints stay hidden unless a test fails.
"""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
TUTORIAL = os.path.join(ROOT, "tutorial")

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


@pytest.fixture
def tutorial():
    """Return a helper that builds paths into the tutorial/ directory."""
    def _path(*parts):
        return os.path.join(TUTORIAL, *parts)
    return _path


@pytest.fixture
def claisen_systems():
    """Parsed g09 ground- and transition-state Systems for the Claisen demo."""
    from pyquiver import quiver
    gs = quiver.System(os.path.join(TUTORIAL, "gaussian", "claisen_gs.out"),
                       style="g09")
    ts = quiver.System(os.path.join(TUTORIAL, "gaussian", "claisen_ts.out"),
                       style="g09")
    return gs, ts
