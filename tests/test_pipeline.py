"""End-to-end test of the titration pipeline on synthetic FCS data.

Generates real FCS files (scatter + GFP) with a known ground-truth titer, writes
the two gate JSONs that ``gate-cells`` would draw interactively, runs the
``analyze`` runner exactly as the CLI does, and asserts the recovered titer,
in-range well count, gating effect, and outputs. Everything is seeded, so it is
deterministic.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import flowio
import numpy as np
import pandas as pd
import pytest

CELLS_SEEDED = 100_000
TRUE_TITER = 3e8  # TU/mL, within the 1e7-1e9 design range
CHANNELS = ["FSC-A", "SSC-A", "FSC-H", "BL1-A"]
DILUTIONS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]  # prep-equivalent uL
N_EVENTS = 20_000


def _make_events(virus_uL: float, rng: np.random.Generator) -> np.ndarray:
    """Scatter with debris + doublets, and a GFP split set by the well's MOI."""
    n = N_EVENTS
    moi = TRUE_TITER * (virus_uL * 1e-3) / CELLS_SEEDED
    frac = 1.0 - np.exp(-moi)
    u = rng.random(n)
    debris = u < 0.12
    doublet = (~debris) & (rng.random(n) < 0.13)
    single = (~debris) & (~doublet)

    fsc_a = np.empty(n); ssc_a = np.empty(n); fsc_h = np.empty(n)
    m = debris; k = int(m.sum())
    fsc_a[m] = rng.normal(6e4, 2e4, k); ssc_a[m] = rng.normal(6e4, 3e4, k)
    fsc_h[m] = fsc_a[m] * 0.9 + rng.normal(0, 1e4, k)
    m = single; k = int(m.sum())
    fa = rng.normal(3e5, 5e4, k)
    fsc_a[m] = fa; ssc_a[m] = rng.normal(2e5, 5e4, k)
    fsc_h[m] = fa * 0.9 + rng.normal(0, 1.5e4, k)
    m = doublet; k = int(m.sum())
    fsc_a[m] = rng.normal(6e5, 7e4, k); ssc_a[m] = rng.normal(2.2e5, 5e4, k)
    fsc_h[m] = rng.normal(2.7e5, 2e4, k)

    fsc_a = np.clip(fsc_a, 1, None); ssc_a = np.clip(ssc_a, 1, None)
    fsc_h = np.clip(fsc_h, 1, None)
    pos = rng.random(n) < frac  # GFP positivity independent of scatter
    gfp = np.where(pos, rng.lognormal(np.log(1e4), 0.40, n),
                   rng.lognormal(np.log(1e2), 0.42, n))
    return np.column_stack([fsc_a, ssc_a, fsc_h, gfp]).astype(np.float32)


def _write_fcs(path: Path, events: np.ndarray) -> None:
    with open(path, "wb") as fh:
        flowio.create_fcs(fh, events.flatten().tolist(), CHANNELS)


def _save_gate(path: Path, x_ch: str, y_ch: str, verts) -> None:
    path.write_text(json.dumps(
        {"name": path.stem.replace("_gate", ""), "x_channel": x_ch, "y_channel": y_ch,
         "vertices": [list(map(float, v)) for v in verts]}))


CONFIG = """\
sheet = "samples.csv"
data  = "data"
out   = "results"
channels = ["BL1-A"]
control = "uninfected"
dose = "virus_uL"

[[analysis]]
kind = "mfi_pct"

[[analysis]]
kind = "titer"
cells_seeded = 100000
volume = "virus_uL"

[[analysis]]
kind = "histograms"
group = "virus_uL"

[[analysis]]
kind = "dose_response"
dose = "virus_uL"
"""


@pytest.fixture(scope="module")
def results_dir(tmp_path_factory) -> Path:
    """Build a synthetic experiment, run the analyze runner, return results/."""
    base = tmp_path_factory.mktemp("titration")
    data = base / "data"; data.mkdir()
    results = base / "results"; results.mkdir()
    rng = np.random.default_rng(20)

    rows = [("uninfected.fcs", "uninfected", 0.0)]
    rows += [(f"v{v:g}.fcs", f"v{v:g}", float(v)) for v in DILUTIONS]
    with open(base / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["filename", "sample", "virus_uL"])
        for fname, label, v in rows:
            _write_fcs(data / fname, _make_events(v, rng))
            w.writerow([fname, label, v])

    _save_gate(results / "cells_gate.json", "FSC-A", "SSC-A",
               [(1.5e5, 5e4), (7.5e5, 5e4), (7.5e5, 4e5), (1.5e5, 4e5)])
    _save_gate(results / "singlets_gate.json", "FSC-A", "FSC-H",
               [(1.5e5, 1.1e5), (5e5, 3.75e5), (5e5, 5.25e5), (1.5e5, 1.6e5)])

    config = base / "titration.toml"
    config.write_text(CONFIG)

    from flowsmith.runner import main
    main([str(config)])
    return results


def test_outputs_written(results_dir):
    for name in ("fluor_stats.csv", "titer.csv", "titer.png",
                 "hist_BL1-A.png", "dose_response.png"):
        assert (results_dir / name).exists(), f"missing output {name}"


def test_gating_removed_debris_and_doublets(results_dir):
    # 12% debris + ~11% doublets removed -> ~76% of N_EVENTS retained as singlets.
    n = pd.read_csv(results_dir / "fluor_stats.csv")["n"]
    assert (n < 0.90 * N_EVENTS).all(), "gating should remove debris/doublets"
    assert (n > 0.60 * N_EVENTS).all(), "gating removed too many events"


def test_recovered_titer_close_to_truth(results_dir):
    titer = pd.read_csv(results_dir / "titer.csv")
    in_range = titer[titer["in_range"]]
    mean = in_range["titer_TU_per_mL"].mean()
    # Poisson-corrected mean should land within ~20% of the ground truth.
    assert abs(mean / TRUE_TITER - 1) < 0.20, f"recovered {mean:.2e} vs {TRUE_TITER:.2e}"


def test_at_least_three_in_range_wells(results_dir):
    titer = pd.read_csv(results_dir / "titer.csv")
    assert int(titer["in_range"].sum()) >= 3


def test_mock_control_is_low_and_saturated_wells_are_nan(results_dir):
    titer = pd.read_csv(results_dir / "titer.csv").set_index("virus_uL")
    assert titer.loc[0.0, "pct_pos_BL1-A"] < 3.0  # mock ~ threshold percentile
    # 100% positive wells cannot give a Poisson titer.
    assert pd.isna(titer.loc[10.0, "titer_TU_per_mL"])
