# *PyQuiver*

[![CI](https://github.com/ekwan/PyQuiver/actions/workflows/ci.yml/badge.svg)](https://github.com/ekwan/PyQuiver/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pyquiver-kie.svg)](https://pypi.org/project/pyquiver-kie/)

*A user-friendly program for calculating isotope effects.*

*PyQuiver* is an open-source Python program for calculating kinetic isotope effects (KIEs) and equilibrium isotope effects (EIEs) using harmonic frequencies and the Bigeleisen-Mayer equation.  *PyQuiver* requires Cartesian Hessian matrices, which can be calculated using any electronic structure program.  

### Features

* automatically read frequencies from [`Gaussian`](https://gaussian.com) (g09/g16) and [`ORCA`](https://orcaforum.kofo.mpg.de/app.php/portal) output files
* excellent performance for larger systems, with optional multi-threading
* custom isotopic substitutions and arbitrary temperature
* tunnelling corrections: Wigner, Bell inverted parabola, and Skodje-Truhlar
* command line **or** a clean Python API with structured (pandas-ready) results

### Installation

*PyQuiver* requires Python 3.9+ and [`numpy`](https://numpy.org).  [`pandas`](https://pandas.pydata.org) is optional (only for `to_dataframe()`).

Install from PyPI:

```bash
pip install pyquiver-kie            # core (numpy only)
pip install "pyquiver-kie[pandas]"  # adds DataFrame output
```

The distribution is named **`pyquiver-kie`** on PyPI, but you still import it as `pyquiver`.

Or install from a checkout:

```bash
git clone https://github.com/ekwan/PyQuiver.git
cd PyQuiver
pip install -e .
```

*PyQuiver* has been tested on PC, Mac, and Linux platforms.

### Command line (deprecated)

The `pyquiver` command still works but is **deprecated** in favour of the Python API below:

```bash
pyquiver demo.config gs.out ts.out            # gaussian by default; prints a deprecation notice
```

Add `-v`/`-vv` for logging, `-s` for the style (`gaussian`/`g16`/`g09`, `orca`, `pyquiver`), `-j` for threads.

### Python API

Build the configuration in Python with `Config.from_dict` (or load a `.config`
file with `Config(path)`), then run it on a ground- and transition-state file:

```python
from pyquiver import Config, KIE_Calculation

config = Config.from_dict(
    isotopologues={"C1": [(1, 1, "13C")], "H/D": [(7, 7, "2D"), (8, 8, "2D")]},
    temperature=393, scaling=0.9614, imag_threshold=50,
)
calc = KIE_Calculation(config, "gs.out", "ts.out", style="gaussian")

print(calc)                       # human-readable table
calc.to_dict()                    # {name: {uncorrected, wigner, inverted_parabola}}
calc.to_csv("kies.csv")           # CSV text / file
calc.to_dataframe()               # pandas DataFrame (needs the pandas extra)
calc.results["C1"].wigner         # named access to each result
```

`System` objects hold a parsed geometry and Hessian, so you can read each file
once and reuse it, which is handy for a parameter scan:

```python
from pyquiver import System

gs = System("gs.out", style="gaussian")   # parse once
ts = System("ts.out", style="gaussian")
for T in (298, 350, 393):
    cfg = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                           temperature=T, scaling=0.9614, imag_threshold=50)
    print(T, KIE_Calculation(cfg, gs, ts).results["C1"].uncorrected)
```

`KIE_Calculation(..., n_jobs=N)` parallelizes the per-isotopologue work across
`N` threads (the heavy `eigvalsh` step releases the GIL). The default is serial;
it pays off for large systems and many isotopologues.

The uncorrected, Wigner, and Bell corrections are reported automatically. For
the Skodje-Truhlar correction, supply the reactant/product/TS single-point
energies (it needs the barrier height):

```python
# energies in hartree by default (or unit="J")
calc.skodje_truhlar(reactant_energy=E_sm, product_energy=E_pr, ts_energy=E_ts)
# -> {isotopologue: corrected_KIE}
```

To run one configuration over many structures, use `batch` with a
`{label: (gs, ts)}` dictionary and get back one table:

```python
from pyquiver import batch

results = batch("demo.config", {
    "b3lyp": ("b3lyp_gs.out", "b3lyp_ts.out"),
    "m06":   ("m06_gs.out",   "m06_ts.out"),
})
results.to_dataframe()        # label, name, uncorrected, wigner, inverted_parabola
results["b3lyp"]              # the KIE_Calculation for that pair
```

You normally build that dictionary from the files on disk, e.g.
`{f.removesuffix("_gs.out"): (f, f.replace("_gs", "_ts")) for f in glob.glob("*_gs.out")}`.
Pass `energies={label: (reactant, ts, product)}` to add a Skodje-Truhlar column.

PyQuiver logs through the standard `logging` module (logger name `pyquiver`) and
raises exceptions on bad input; it never prints on import or calls `sys.exit`,
so it is safe to use inside scripts and notebooks.

## Tutorial

To learn how to use *PyQuiver*, please look at the [tutorial](tutorial/TUTORIAL.md).

## Further Reading

* [Technical Details](DETAILS.md) (compatibility with QUIVER, how to invoke PyQuiver, `.config` files, input files, the PyQuiver Standard format, miscellaneous notes)
* [Batch processing](src/AUTOQUIVER.md) (how to run PyQuiver on many related files)
* [Snipping Utility](scripts/SNIP.md) (how to "snip" out only the relevant parts of a Gaussian output file for publication or archiving purposes)

## Authors

*PyQuiver* was written by Thayer Anderson and Eugene Kwan.  Please email `ekwan16@gmail.com` with any questions.  We will gladly try to help you.

We thank Gregor Giesen for writing the code to read ORCA output.  We thank Corin Wagen for miscellaneous optimization.  We also thank Christoph Riplinger and Giovanni Bistoni for valuable discussions.

## How to Cite

Anderson, T.L.; Kwan, E.E.  *PyQuiver*, `www.github.com/ekwan/PyQuiver` (`pip install pyquiver-kie`)

## License
   
This project is licensed under the Apache License, Version 2.0. See `LICENSE.txt` for full terms and conditions.
   
*Copyright 2020 by Thayer L. Anderson and Eugene E. Kwan*
