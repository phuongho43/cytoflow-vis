"""Tests for the automatic gating heuristics."""
from __future__ import annotations

import numpy as np

from flowsmith.gating import auto_ellipse_gate, auto_singlet_gate, points_in_polygon
import pandas as pd


def _inside(verts, x, y, x_ch="x", y_ch="y"):
    df = pd.DataFrame({x_ch: x, y_ch: y})
    return points_in_polygon(df, x_ch, y_ch, verts)


def test_ellipse_gate_keeps_cells_drops_debris():
    rng = np.random.default_rng(0)
    cells = rng.multivariate_normal([3e5, 2e5], [[5e4**2, 0], [0, 5e4**2]], 8000)
    debris = rng.multivariate_normal([6e4, 6e4], [[2e4**2, 0], [0, 2e4**2]], 2000)
    x = np.concatenate([cells[:, 0], debris[:, 0]])
    y = np.concatenate([cells[:, 1], debris[:, 1]])

    verts = auto_ellipse_gate(x, y)
    inside = _inside(verts, x, y)
    cell_mask = np.arange(len(x)) < len(cells)

    assert inside[cell_mask].mean() > 0.9, "should keep most cells"
    assert inside[~cell_mask].mean() < 0.02, "should drop the debris cluster"


def test_singlet_gate_keeps_singlets_drops_doublets():
    rng = np.random.default_rng(1)
    fsc_a_s = rng.normal(3e5, 5e4, 8000)
    singlets = np.column_stack([fsc_a_s, fsc_a_s * 0.9 + rng.normal(0, 1.5e4, 8000)])
    # doublets: ~2x area, single-cell height -> low FSC-H/FSC-A ratio
    doublets = np.column_stack([rng.normal(6e5, 7e4, 2000), rng.normal(2.7e5, 2e4, 2000)])
    x = np.concatenate([singlets[:, 0], doublets[:, 0]])
    y = np.concatenate([singlets[:, 1], doublets[:, 1]])

    verts = auto_singlet_gate(x, y)
    inside = _inside(verts, x, y)
    singlet_mask = np.arange(len(x)) < len(singlets)

    assert inside[singlet_mask].mean() > 0.9, "should keep most singlets"
    assert inside[~singlet_mask].mean() < 0.05, "should drop doublets (off the diagonal)"
