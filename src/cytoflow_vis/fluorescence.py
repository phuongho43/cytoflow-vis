"""Fluorescence analysis of a gated (singlet) population.

Works on the *populations* produced by the gating pipeline (a list of dicts
with ``sample_id``, ``conditions`` and an ``events`` DataFrame). Fluorescence
channels are displayed with FlowKit's logicle transform, which maps raw values
into ~[0, 1] display units while handling negatives and the compressed low end.

Provides: per-sample stats (MFI + % positive), 1D histograms grouped by
condition, 2D channel-vs-channel density, and dose-response plots.
"""

from __future__ import annotations

import flowkit as fk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cytoflow_vis.plotting import density_plot

DEFAULT_FLUOR_CHANNELS = ["BL1-A", "RL1-A"]

# Logicle display range (a touch below 0 to show the negative population).
_LOGICLE_LO, _LOGICLE_HI = -0.15, 1.05


def make_logicle(
    param_t: float = 262144.0,
    param_w: float = 0.5,
    param_m: float = 4.5,
    param_a: float = 0.0,
) -> fk.transforms.LogicleTransform:
    """Standard logicle transform. ``param_t`` is the top of the raw scale."""
    return fk.transforms.LogicleTransform(
        param_t=param_t, param_w=param_w, param_m=param_m, param_a=param_a
    )


def _subsample(values: np.ndarray, n: int, rng) -> np.ndarray:
    if n is None or len(values) <= n:
        return values
    return values[rng.choice(len(values), size=n, replace=False)]


def compute_stats(
    populations: list[dict],
    channels: list[str],
    xform,
    control_id: str | None = None,
    positive_percentile: float = 99.0,
) -> tuple[pd.DataFrame, dict]:
    """Per-sample fluorescence stats.

    Returns ``(stats_df, thresholds)``. MFI is the median of the *raw* values
    (standard reporting). If ``control_id`` is given, a positive threshold is
    set per channel at the given percentile of the control's transformed
    values, and ``pct_pos_<channel>`` is the fraction of events above it.
    """
    thresholds: dict[str, float] = {}
    if control_id is not None:
        control = next((p for p in populations if p["sample_id"] == control_id), None)
        if control is None:
            raise ValueError(f"control sample {control_id!r} not found among populations")
        for ch in channels:
            xf = xform.apply(control["events"][ch].to_numpy(dtype=float))
            thresholds[ch] = float(np.percentile(xf, positive_percentile))

    rows = []
    for p in populations:
        row = {"sample_id": p["sample_id"], **p["conditions"], "n": int(len(p["events"]))}
        for ch in channels:
            raw = p["events"][ch].to_numpy(dtype=float)
            row[f"MFI_{ch}"] = round(float(np.median(raw)), 2) if len(raw) else float("nan")
            if ch in thresholds:
                xf = xform.apply(raw)
                row[f"pct_pos_{ch}"] = (
                    round(100.0 * float(np.mean(xf > thresholds[ch])), 2) if len(xf) else float("nan")
                )
        rows.append(row)
    return pd.DataFrame(rows), thresholds


def plot_histograms(
    populations: list[dict],
    channel: str,
    xform,
    ax: plt.Axes | None = None,
    group_col: str | None = None,
    bins: int = 200,
    per_sample: int = 20000,
    threshold: float | None = None,
    seed: int = 0,
) -> plt.Axes:
    """Overlaid logicle histograms of ``channel``, one line per sample."""
    rng = np.random.default_rng(seed)
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))

    def sort_key(p):
        val = p["conditions"].get(group_col) if group_col else p["sample_id"]
        return (val is None, val)

    pops = sorted(populations, key=sort_key)
    cmap = plt.get_cmap("viridis")
    for i, p in enumerate(pops):
        raw = _subsample(p["events"][channel].to_numpy(dtype=float), per_sample, rng)
        xf = xform.apply(raw)
        color = cmap(i / max(1, len(pops) - 1))
        label = p["sample_id"]
        if group_col:
            label = f"{p['sample_id']} ({group_col}={p['conditions'].get(group_col)})"
        ax.hist(
            xf, bins=bins, range=(_LOGICLE_LO, _LOGICLE_HI), histtype="step",
            density=True, color=color, lw=1.5, label=label,
        )
    if threshold is not None:
        ax.axvline(threshold, color="k", ls="--", lw=1, label="+ threshold")
    ax.set_xlabel(f"{channel} (logicle)")
    ax.set_ylabel("density")
    ax.set_title(f"{channel} distribution by sample")
    ax.legend(fontsize=8)
    return ax


def plot_2d_density(
    populations: list[dict],
    x_channel: str,
    y_channel: str,
    xform,
    per_sample: int = 20000,
    ncols: int = 3,
    bins: int = 150,
    seed: int = 0,
) -> plt.Figure:
    """Faceted logicle 2D density (one panel per sample)."""
    rng = np.random.default_rng(seed)
    n = len(populations)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
    flat = axes.flatten()
    for ax, p in zip(flat, populations):
        x = xform.apply(_subsample(p["events"][x_channel].to_numpy(dtype=float), per_sample, rng))
        y = xform.apply(_subsample(p["events"][y_channel].to_numpy(dtype=float), per_sample, rng))
        density_plot(
            x, y, ax=ax, xlabel=f"{x_channel} (logicle)", ylabel=f"{y_channel} (logicle)",
            title=p["sample_id"], bins=bins, clip_percentile=None,
        )
        ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
        ax.set_ylim(_LOGICLE_LO, _LOGICLE_HI)
    for ax in flat[n:]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def plot_dose_response(
    ax: plt.Axes,
    stats_df: pd.DataFrame,
    dose_col: str,
    y_col: str,
    group_col: str | None = None,
    logx: bool = True,
) -> plt.Axes:
    """Plot ``y_col`` vs ``dose_col``, one line per ``group_col`` value."""
    if group_col and group_col in stats_df and stats_df[group_col].nunique() > 1:
        for gval, sub in stats_df.groupby(group_col):
            sub = sub.sort_values(dose_col)
            ax.plot(sub[dose_col], sub[y_col], marker="o", label=f"{group_col}={gval}")
        ax.legend(fontsize=8)
    else:
        sub = stats_df.sort_values(dose_col)
        ax.plot(sub[dose_col], sub[y_col], marker="o")
    if logx:
        ax.set_xscale("symlog")  # symlog tolerates a zero dose
    ax.set_ylim(bottom=0)  # anchor at 0 so a flat channel reads as flat
    ax.set_xlabel(dose_col)
    ax.set_ylabel(y_col)
    return ax
