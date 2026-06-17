"""Deprecated. PyQuiver now uses the standard ``logging`` module instead of a
global DEBUG level.

To see diagnostics, configure logging in your own code, e.g.::

    import logging
    logging.basicConfig(level=logging.INFO)   # or DEBUG for more detail

``DEBUG`` is kept only so that old code importing it does not break; it no
longer controls any output.
"""

DEBUG = 0
