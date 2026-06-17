# Batch processing

To run one configuration over many ground-state/transition-state pairs, use
`pyquiver.batch` with a `{label: (gs, ts)}` dictionary. File discovery is left
to ordinary Python, so no filename convention is assumed; `batch` just maps the
configuration over the pairs you give it and collects the results into a table.

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

## Building the pairs

You usually build the dictionary from the files on disk. Glob the ground-state
files, derive a label from each filename, and pair it with the matching
transition state:

```python
import glob
from pyquiver import batch

label = lambda path: path.removesuffix("_gs.out")
pairs = {label(gs): (gs, gs.replace("_gs", "_ts")) for gs in sorted(glob.glob("*_gs.out"))}
results = batch("demo.config", pairs)
```

Adjust the `label`/`replace` rule to match your own naming.

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
