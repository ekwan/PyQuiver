# Changelog

All notable changes to *PyQuiver* are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [1.1.0] - 2026-06-17

First PyPI release of *PyQuiver* as an installable package
(`pip install pyquiver-kie`; imported as `pyquiver`).

### Added
- Installable `pyquiver` package with a clean, scriptable Python API:
  `Config.from_dict`, `KIE_Calculation`, `System`, `Isotopologue`, and `batch`.
- Structured, pandas-ready results (`to_dict`, `to_csv`, `to_dataframe`, and
  named access such as `calc.results["C1"].wigner`).
- Skodje-Truhlar tunnelling correction (`calc.skodje_truhlar(...)`), alongside
  the existing Wigner and Bell (infinite parabola) corrections.
- Support for custom and unusual isotope masses, including bare numeric masses.
- Optional multi-threading via `n_jobs` for large systems and many
  isotopologues (the heavy `eigvalsh` step releases the GIL).
- A warning when predicting primary H/D KIEs (unreliable due to tunnelling) and
  when the Skodje-Truhlar correction is requested below its crossover
  temperature.
- ORCA `.hess` parsing in addition to Gaussian (g09/g16) and the native `.qin`
  format.
- Full test suite across Python 3.9-3.13, GitHub Actions CI, a publish workflow
  using PyPI Trusted Publishing, and a rewritten tutorial notebook.

### Changed
- Reorganised the source tree into an installable package layout.
- Renamed the Bell correction terminology to "infinite parabola" throughout the
  code, results, and documentation.

### Deprecated
- The `pyquiver` command line is deprecated in favour of the Python API. It
  still works but prints a deprecation notice.

[1.1.0]: https://github.com/ekwan/PyQuiver/releases/tag/v1.1.0
