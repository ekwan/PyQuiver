"""Command-line entry point for PyQuiver (the ``pyquiver`` console script).

The command line is deprecated; use the Python API (pyquiver.KIE_Calculation)
instead. Batch processing over many files is Python-only; see pyquiver.batch.
"""

import sys
import argparse
import logging
import warnings

from .kie import KIE_Calculation


def _configure_logging(verbosity):
    """Map a -v count to a logging level and install a simple handler.

    0 -> WARNING (quiet), 1 -> INFO, 2+ -> DEBUG.
    """
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity and verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(message)s")


def main(argv=None):
    """Run a single KIE/EIE calculation (the ``pyquiver`` command)."""
    parser = argparse.ArgumentParser(
        description="PyQuiver calculates KIEs and EIEs based on a ground "
                    "and transition state file.")
    parser.add_argument('-v', '--verbose', dest="debug", action='count',
                        help='print debug information (repeatable)')
    parser.add_argument('-s', '--style', dest="style", default='gaussian',
                        help='style of input files (gaussian, orca, or pyquiver)')
    parser.add_argument('-j', '--jobs', dest="jobs", type=int, default=1,
                        help='number of worker threads for the isotopologue '
                             'computations (default 1; -1 uses all cores)')
    parser.add_argument('config', help='configuration file path')
    parser.add_argument('gs', help='ground state file path')
    parser.add_argument('ts', help='transition state file path')

    args = parser.parse_args(argv)
    _configure_logging(args.debug or 0)

    _msg = ("the 'pyquiver' command line is deprecated and will be removed in a "
            "future release; use the Python API (pyquiver.KIE_Calculation) instead.")
    warnings.warn(_msg, DeprecationWarning, stacklevel=2)
    print("warning: " + _msg, file=sys.stderr)

    calc = KIE_Calculation(args.config, args.gs, args.ts, style=args.style,
                           n_jobs=args.jobs)
    print(calc)
    return calc
