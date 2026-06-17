"""Run a KIE/EIE calculation over many ground-state/transition-state pairs.

File discovery and pairing are left to ordinary Python (glob, pathlib, ...);
``batch`` just maps one configuration over the pairs you give it and collects
the results into a table. No filename convention is assumed::

    import glob
    from pyquiver import batch

    gs = sorted(glob.glob("*_reactant.out"))
    ts = sorted(glob.glob("*_ts.out"))
    results = batch("demo.config", zip(gs, ts))   # labels inferred from filenames
    results.to_dataframe()

Pairs may be given as a mapping ``{label: (gs, ts)}``, an iterable of
``(label, gs, ts)`` triples, or an iterable of ``(gs, ts)`` pairs (in which
case the label is the ground-state filename without its extension).
"""

import os
from collections import OrderedDict

from .config import Config
from .kie import KIE_Calculation


def _normalize_pairs(pairs):
    """Yield (label, gs, ts) from any of the supported pair formats."""
    items = pairs.items() if hasattr(pairs, "items") else pairs
    for item in items:
        if hasattr(pairs, "items"):
            label, (gs, ts) = item
        elif len(item) == 3:
            label, gs, ts = item
        elif len(item) == 2:
            gs, ts = item
            label = os.path.splitext(os.path.basename(str(gs)))[0]
        else:
            raise ValueError("each pair must be (gs, ts) or (label, gs, ts)")
        yield label, gs, ts


class BatchResults(object):
    """Results of a :func:`batch` run: one KIE_Calculation per label."""

    def __init__(self, calcs, skodje_truhlar=None):
        self._calcs = OrderedDict(calcs)   # label -> KIE_Calculation
        self._st = skodje_truhlar          # label -> {isotopologue: corrected} or None

    def __getitem__(self, label):
        return self._calcs[label]

    def __iter__(self):
        return iter(self._calcs)

    def __len__(self):
        return len(self._calcs)

    def items(self):
        return self._calcs.items()

    def to_records(self):
        """Tidy long-form rows: one dict per (label, isotopologue)."""
        rows = []
        for label, calc in self._calcs.items():
            for result in calc.results:
                row = {"label": label}
                row.update(result._asdict())
                if self._st is not None:
                    row["skodje_truhlar"] = self._st[label].get(result.name)
                rows.append(row)
        return rows

    def to_dataframe(self):
        """Return the results as a pandas DataFrame (optional pandas extra)."""
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("to_dataframe requires pandas; install with "
                              "`pip install pyquiver-kie[pandas]`")
        return pd.DataFrame(self.to_records())

    def to_csv(self, path=None, float_format="%.4f"):
        """Render the results as CSV text; write to ``path`` if given."""
        records = self.to_records()
        if not records:
            return ""
        columns = list(records[0].keys())

        def cell(v):
            return float_format % v if isinstance(v, float) else str(v)

        lines = [",".join(columns)]
        lines += [",".join(cell(rec[c]) for c in columns) for rec in records]
        text = "\n".join(lines) + "\n"
        if path is not None:
            with open(path, "w") as f:
                f.write(text)
        return text


def batch(config, pairs, style="gaussian", n_jobs=1, energies=None):
    """Run a KIE calculation for each ground-state/transition-state pair.

    ``config`` is a path to a .config file or a :class:`~pyquiver.Config`.
    ``pairs`` is a mapping ``{label: (gs, ts)}`` or an iterable of ``(gs, ts)``
    or ``(label, gs, ts)`` (see the module docstring). Returns
    :class:`BatchResults`.

    If ``energies`` is given, a Skodje-Truhlar correction is computed for every
    pair and included in the results. It maps each label to the three
    single-point electronic energies the barrier needs, in hartree, ordered
    along the reaction coordinate:
    ``{label: (reactant_energy, ts_energy, product_energy)}``. (If your reaction
    is effectively a single well, pass the reactant energy for the product too.)
    """
    if isinstance(config, str):
        config = Config(config)   # parse once, reuse for every pair

    calcs = OrderedDict()
    st = OrderedDict() if energies is not None else None
    for label, gs, ts in _normalize_pairs(pairs):
        calc = KIE_Calculation(config, gs, ts, style=style, n_jobs=n_jobs)
        calcs[label] = calc
        if energies is not None:
            reactant, ts_energy, product = energies[label]
            st[label] = calc.skodje_truhlar(reactant, product, ts_energy)
    return BatchResults(calcs, st)
