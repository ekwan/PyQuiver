"""Structured result objects for KIE/EIE calculations.

These give named access to results instead of positional indexing into
``KIE.value``, and convert a whole calculation to a dict, CSV, or pandas
DataFrame. pandas is optional: ``to_dict``/``to_csv`` are dependency-free,
``to_dataframe`` requires ``pandas`` (``pip install pyquiver-kie[pandas]``).
"""

from collections import namedtuple

# named view of a kinetic isotope effect (the three reported columns)
KIEResult = namedtuple("KIEResult",
                       ["name", "uncorrected", "wigner", "inverted_parabola"])

# named view of an equilibrium isotope effect
EIEResult = namedtuple("EIEResult", ["name", "value"])


class Results(object):
    """Ordered, tabular view of a KIE_Calculation's isotopologue results."""

    def __init__(self, calc):
        self.is_eie = (calc.eie_flag == 1)
        self.reference_isotopologue = calc.config.reference_isotopologue
        self._items = [k.result for k in calc.KIES.values()]

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)

    def __getitem__(self, name):
        for r in self._items:
            if r.name == name:
                return r
        raise KeyError(name)

    @property
    def columns(self):
        return ["name", "value"] if self.is_eie else \
            ["name", "uncorrected", "wigner", "inverted_parabola"]

    def to_records(self):
        """Return a list of dicts, one per isotopologue."""
        return [r._asdict() for r in self._items]

    def to_dict(self):
        """Map isotopologue name -> value (EIE) or column dict (KIE)."""
        if self.is_eie:
            return {r.name: r.value for r in self._items}
        return {r.name: {"uncorrected": r.uncorrected,
                         "wigner": r.wigner,
                         "inverted_parabola": r.inverted_parabola}
                for r in self._items}

    def to_dataframe(self):
        """Return a pandas DataFrame (requires the optional pandas extra)."""
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("to_dataframe requires pandas; install with "
                              "`pip install pyquiver-kie[pandas]`")
        return pd.DataFrame(self.to_records(), columns=self.columns)

    def to_csv(self, path=None, float_format="%.4f"):
        """Render results as CSV text; write to ``path`` if given, else return."""
        lines = [",".join(self.columns)]
        for r in self._items:
            cells = [r.name] + [float_format % v for v in tuple(r)[1:]]
            lines.append(",".join(cells))
        text = "\n".join(lines) + "\n"
        if path is not None:
            with open(path, "w") as f:
                f.write(text)
        return text
