# Batch processing

To run one configuration over many ground-state/transition-state pairs, use
`pyquiver.batch`. File discovery and pairing are left to ordinary Python, so no
filename convention is assumed; `batch` just maps the configuration over the
pairs you give it and collects the results into a table.

```python
from pyquiver import batch

results = batch("demo.config", {
    "b3lyp": ("b3lyp_gs.out", "b3lyp_ts.out"),
    "m06":   ("m06_gs.out",   "m06_ts.out"),
}, style="gaussian")

results.to_dataframe()          # label, name, uncorrected, wigner, inverted_parabola
results.to_csv("kies.csv")
results["b3lyp"]                # the KIE_Calculation for that pair
```

## Specifying pairs

`pairs` may be given three ways:

* a mapping `{label: (gs, ts)}`
* an iterable of `(label, gs, ts)` triples
* an iterable of `(gs, ts)` pairs — the label is inferred from the ground-state
  filename (without its extension)

So you can discover and pair files however you like, then hand the result to
`batch`. For example, zipping two globs:

```python
import glob
from pyquiver import batch

gs = sorted(glob.glob("*_reactant.out"))
ts = sorted(glob.glob("*_ts.out"))
results = batch("demo.config", zip(gs, ts))
```

## Skodje-Truhlar in batch

Pass `energies` mapping each label to its three single-point electronic
energies (hartree), ordered along the reaction coordinate
`(reactant, ts, product)`, to add a Skodje-Truhlar column:

```python
results = batch(
    "demo.config",
    {"b3lyp": ("b3lyp_gs.out", "b3lyp_ts.out")},
    energies={"b3lyp": (E_reactant, E_ts, E_product)},
)
```

(For distinct products use the real product energy; for a single-well reaction
pass the reactant energy for the product too.)
