"""Input-file parsers and the style registry.

Each parser exposes ``parse(path) -> ParsedSystem``. Styles are looked up
through ``_ALIASES`` so that, e.g., ``gaussian``, ``g16`` and ``g09`` all map
to the Gaussian parser. ``gaussian`` is the preferred name; ``g09`` is kept as
an alias for backward compatibility.
"""

from . import gaussian, orca, native
from ._common import ParsedSystem

__all__ = ["ParsedSystem", "parse", "supported_styles"]

_PARSERS = {
    "gaussian": gaussian.parse,
    "orca": orca.parse,
    "native": native.parse,
}

# user-facing style name -> canonical parser key
_ALIASES = {
    "gaussian": "gaussian",
    "g16": "gaussian",
    "g09": "gaussian",      # deprecated spelling, still accepted
    "orca": "orca",
    "native": "native",
    "pyquiver": "native",   # the .qin format
    "qin": "native",
}


def supported_styles():
    """Return the sorted list of accepted style names."""
    return sorted(_ALIASES)


def parse(path, style):
    """Parse ``path`` using the parser registered for ``style``."""
    key = _ALIASES.get(style)
    if key is None:
        raise ValueError("specified style, {0}, not supported (choose from {1})"
                         .format(style, ", ".join(supported_styles())))
    return _PARSERS[key](path)
