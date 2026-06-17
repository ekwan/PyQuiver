# PyQuiver test suite

Standard-library `unittest` tests (no extra dependencies beyond NumPy, which
PyQuiver already requires). Run from the repository root:

```bash
python3 -m unittest discover -s test
```

or run a single module / verbose:

```bash
python3 -m unittest test.test_regression -v
```

`pytest` also works if installed (`pytest test`).

## What is covered

| File | Scope |
|------|-------|
| `test_regression.py` | End-to-end KIE values for the Claisen tutorial (g09 + ORCA inputs) locked to known-good numbers. The safety net for refactoring. |
| `test_physics.py` | Pure functions: `u`, `wigner`, `bell`, `partition_components`. |
| `test_isotopologue.py` | Mass-weighting, frequency caching, and the reference-isotopologue caching invariant (the reference is diagonalized once, not once per substitution). |
| `test_parsers.py` | `quiver.System` parsing for g09, ORCA, and the native `.qin` format (including a round-trip check against the g09 Hessian). |

`context.py` adds `src/` to `sys.path` and provides path/quiet helpers.
