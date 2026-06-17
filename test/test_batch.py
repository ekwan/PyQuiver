"""Tests for the batch runner (pyquiver.batch)."""

import pytest

from pyquiver import batch, KIE_Calculation
from pyquiver.batch import BatchResults

CONFIG = ("gaussian", "claisen_demo.config")
GS = ("gaussian", "claisen_gs.out")
TS = ("gaussian", "claisen_ts.out")


@pytest.fixture
def files(tutorial):
    return tutorial(*CONFIG), tutorial(*GS), tutorial(*TS)


def test_pairs_as_mapping(files):
    cfg, gs, ts = files
    results = batch(cfg, {"a": (gs, ts), "b": (gs, ts)})
    assert set(results) == {"a", "b"}
    assert len(results) == 2
    assert "C1" in results["a"].KIES        # __getitem__ -> KIE_Calculation
    assert [label for label, _ in results.items()] == ["a", "b"]   # order preserved


def test_bad_pair_value_rejected(files):
    cfg, gs, ts = files
    with pytest.raises(ValueError):
        batch(cfg, {"a": (gs,)})            # value must be a (gs, ts) pair


def test_to_records_and_dataframe(files):
    cfg, gs, ts = files
    results = batch(cfg, {"run": (gs, ts)})
    records = results.to_records()
    assert {r["label"] for r in records} == {"run"}
    assert "uncorrected" in records[0]
    pd = pytest.importorskip("pandas")
    df = results.to_dataframe()
    assert "label" in df.columns and len(df) == len(records)


def test_to_csv(files, tmp_path):
    cfg, gs, ts = files
    results = batch(cfg, {"run": (gs, ts)})
    path = tmp_path / "out.csv"
    text = results.to_csv(str(path))
    assert path.read_text() == text
    assert text.splitlines()[0].startswith("label,")


def test_empty_batch_csv():
    assert BatchResults({}).to_csv() == ""


def test_to_dataframe_without_pandas(files, monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "pandas":
            raise ImportError("no pandas")
        return real_import(name, *a, **k)

    cfg, gs, ts = files
    results = batch(cfg, {"run": (gs, ts)})
    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError):
        results.to_dataframe()


def test_skodje_truhlar_column(files):
    cfg, gs, ts = files
    # energies per label: (reactant, ts, product) in hartree
    results = batch(cfg, {"run": (gs, ts)},
                    energies={"run": (0.0, 0.02, 0.0)})
    by_iso = {r["name"]: r for r in results.to_records()}
    assert "skodje_truhlar" in by_iso["C1"]
    assert by_iso["C1"]["skodje_truhlar"] is not None

    # the batch energy tuple is (reactant, ts, product); confirm it is mapped
    # onto skodje_truhlar(reactant, product, ts) and not silently transposed
    direct = KIE_Calculation(cfg, gs, ts).skodje_truhlar(
        reactant_energy=0.0, product_energy=0.0, ts_energy=0.02)
    assert by_iso["C1"]["skodje_truhlar"] == pytest.approx(direct["C1"])
