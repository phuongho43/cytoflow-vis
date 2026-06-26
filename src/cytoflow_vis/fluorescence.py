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
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator

from cytoflow_vis.plotting import density_plot
from cytoflow_vis.style import CATEGORICAL_PALETTE, INK, sequential_colors

DEFAULT_FLUOR_CHANNELS = ["BL1-A", "RL1-A"]

# Logicle display range (a touch below 0 to show the negative population).
_LOGICLE_LO, _LOGICLE_HI = -0.15, 1.05

# Decades labelled on a logicle (biexponential) axis, the cytometry convention.
_LOGICLE_DECADES = [100, 1000, 10000, 100000]


def logicle_axis(ax, xform, which: str = "x", lo: float = _LOGICLE_LO, hi: float = _LOGICLE_HI) -> None:
    """Relabel a logicle display axis with biexponential decade ticks.

    The events are plotted in logicle *display* units (~[0, 1]); this maps the
    ticks back to raw values (0, 10^2 ... 10^5) at their display positions, with
    minor ticks at the 2-9 subdivisions, so the axis reads like FlowJo/Attune.
    """
    major_raw = np.array([0, *_LOGICLE_DECADES], dtype=float)
    minor_raw = np.array(
        [k * d for d in [10, *_LOGICLE_DECADES] for k in range(2, 10)], dtype=float
    )
    major_pos = xform.apply(major_raw)
    minor_pos = xform.apply(minor_raw)
    labels = ["0", *[f"$10^{int(round(np.log10(d)))}$" for d in _LOGICLE_DECADES]]

    in_range = (major_pos >= lo) & (major_pos <= hi)
    major_pos, labels = major_pos[in_range], [l for l, keep in zip(labels, in_range) if keep]
    minor_pos = minor_pos[(minor_pos >= lo) & (minor_pos <= hi)]

    axis = ax.xaxis if which == "x" else ax.yaxis
    axis.set_major_locator(FixedLocator(major_pos))
    axis.set_major_formatter(FixedFormatter(labels))
    axis.set_minor_locator(FixedLocator(minor_pos))


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


def control_thresholds(
    populations: list[dict],
    channels: list[str],
    xform,
    control_id: str,
    percentile: float = 99.0,
) -> dict[str, float]:
    """Positive threshold per channel = ``percentile`` of the control's
    transformed values. Use an unstained / untreated control so the threshold
    lands in the negative population."""
    control = next((p for p in populations if p["sample_id"] == control_id), None)
    if control is None:
        raise ValueError(f"control sample {control_id!r} not found among populations")
    return {
        ch: float(np.percentile(xform.apply(control["events"][ch].to_numpy(dtype=float)), percentile))
        for ch in channels
    }


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
    thresholds = (
        control_thresholds(populations, channels, xform, control_id, positive_percentile)
        if control_id is not None
        else {}
    )

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


def _smooth(h: np.ndarray, sigma_bins: float) -> np.ndarray:
    """Gaussian smooth of a 1D histogram (``sigma_bins`` in bin units).

    Edge-corrected: the kernel is renormalised by its own truncated weight near
    the array ends, so the first/last bins are not attenuated. ``sigma_bins <= 0``
    disables smoothing.
    """
    if not sigma_bins or sigma_bins <= 0:
        return h
    radius = int(np.ceil(3 * sigma_bins))
    x = np.arange(-radius, radius + 1)
    k = np.exp(-0.5 * (x / sigma_bins) ** 2)
    k /= k.sum()
    num = np.convolve(h, k, mode="same")
    denom = np.convolve(np.ones_like(h), k, mode="same")  # edge correction
    return num / denom


def n_groups(populations: list[dict], group_col: str | None) -> int:
    """Number of ridges a histogram of these populations will draw."""
    if group_col:
        return len({p["conditions"].get(group_col) for p in populations})
    return len(populations)


def plot_histograms(
    populations: list[dict],
    channel: str,
    xform,
    ax: plt.Axes | None = None,
    group_col: str | None = None,
    group_label: str | None = None,
    colors: list | None = None,
    bins: int = 120,
    per_sample: int = 20000,
    threshold: float | None = None,
    overlap: float = 1.7,
    smooth: float = 1.5,
    seed: int = 0,
) -> plt.Axes:
    """Ridgeline of ``channel`` over a biexponential logicle x-axis.

    One offset, filled ridge per ``group_col`` value (replicates pooled), or per
    sample when no group is given, coloured from the signature ``PALETTE`` and
    outlined in ink. Each ridge's tick is its bare condition value (number or
    string); ``group_label`` names the y-axis (e.g. ``"Dose (mM)"``, defaulting
    to ``group_col``). When ``threshold`` is given, the positive gate is drawn
    and each ridge is annotated with its % positive. Expects the ``rc()`` context.
    """
    rng = np.random.default_rng(seed)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 8))

    # Collect events per ridge: pool replicates sharing a group value.
    def key_of(p):
        return p["conditions"].get(group_col) if group_col else p["sample_id"]

    grouped: dict = {}
    for p in populations:
        grouped.setdefault(key_of(p), []).append(p)
    order = sorted(grouped, key=lambda v: (v is None, v))
    n = len(order)
    step = 1.0
    baselines = []

    # Colour by data type: ordered/numeric conditions get the sequential ramp,
    # categorical (string) conditions get the categorical palette. ``colors``
    # overrides this when given.
    if colors is None:
        numeric = all(
            isinstance(v, (int, float, np.number)) and not isinstance(v, bool) for v in order
        )
        colors = (
            sequential_colors(n)
            if numeric
            else [CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)] for i in range(n)]
        )

    # Draw top row first; later (lower) rows get higher zorder so they occlude.
    for rank, key in enumerate(order):
        baseline = (n - 1 - rank) * step
        baselines.append(baseline)
        color = colors[rank]
        raw = np.concatenate([g["events"][channel].to_numpy(dtype=float) for g in grouped[key]])
        xf = xform.apply(_subsample(raw, per_sample, rng))

        h, edges = np.histogram(xf, bins=bins, range=(_LOGICLE_LO, _LOGICLE_HI), density=True)
        centers = (edges[:-1] + edges[1:]) / 2.0
        h = _smooth(h, smooth)
        h = h / h.max() if h.max() > 0 else h  # per-ridge height normalisation
        curve = baseline + h * overlap
        z = 2 + rank

        ax.fill_between(centers, baseline, curve, facecolor=color, alpha=0.85,
                        edgecolor="none", zorder=z)
        ax.plot(centers, curve, color=INK, lw=1.8, zorder=z,
                path_effects=[pe.Stroke(linewidth=3.2, foreground="white"), pe.Normal()])
        ax.hlines(baseline, _LOGICLE_LO, _LOGICLE_HI, color=INK, lw=1.0, zorder=z)

        if threshold is not None:
            pct = 100.0 * float(np.mean(xf > threshold))
            ax.text(_LOGICLE_HI - 0.01, baseline + 0.12 * overlap, f"{pct:.0f}%",
                    ha="right", va="bottom", color=color, fontweight="bold",
                    fontsize=plt.rcParams["xtick.labelsize"], zorder=n + 5)

    if threshold is not None:
        ax.axvline(threshold, color=INK, ls=(0, (4, 3)), lw=2.5, zorder=n + 4)

    logicle_axis(ax, xform, which="x")
    ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
    ax.set_ylim(-0.25, (n - 1) * step + overlap + 0.35)
    ax.set_yticks(baselines)
    ax.set_yticklabels([str(k) for k in order])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(channel)
    ax.set_ylabel(group_label or group_col or "Sample")
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


# --- Quadrant analysis (two thresholds -> four populations) ------------------
# Quadrant keys, with their corner on an (x, y) plot:
#   dn  = x-/y-  (lower-left)    x_pos = x+/y-  (lower-right)
#   dp  = x+/y+  (upper-right)   y_pos = x-/y+  (upper-left)
QUADRANT_KEYS = ("dn", "x_pos", "dp", "y_pos")


def _resolve_labels(labels: dict | None) -> dict:
    return {**{k: k for k in QUADRANT_KEYS}, **(labels or {})}


def quadrant_stats(
    populations: list[dict],
    x_channel: str,
    y_channel: str,
    xform,
    x_threshold: float,
    y_threshold: float,
    labels: dict | None = None,
) -> pd.DataFrame:
    """Per-sample percentage of events in each of the four quadrants.

    ``labels`` optionally renames the quadrant keys (e.g. for an apoptosis
    panel: ``{"dn": "live", "x_pos": "early", "dp": "late", "y_pos": "necrotic"}``).
    """
    labels = _resolve_labels(labels)
    rows = []
    for p in populations:
        xf = xform.apply(p["events"][x_channel].to_numpy(dtype=float))
        yf = xform.apply(p["events"][y_channel].to_numpy(dtype=float))
        xp, yp = xf > x_threshold, yf > y_threshold
        n = len(xf)
        frac = {
            "dn": np.mean(~xp & ~yp) if n else float("nan"),
            "x_pos": np.mean(xp & ~yp) if n else float("nan"),
            "dp": np.mean(xp & yp) if n else float("nan"),
            "y_pos": np.mean(~xp & yp) if n else float("nan"),
        }
        row = {"sample_id": p["sample_id"], **p["conditions"], "n": int(n)}
        for k in QUADRANT_KEYS:
            row[f"pct_{labels[k]}"] = round(100.0 * float(frac[k]), 2)
        rows.append(row)
    return pd.DataFrame(rows)


def plot_quadrants(
    populations: list[dict],
    x_channel: str,
    y_channel: str,
    xform,
    x_threshold: float,
    y_threshold: float,
    quad_df: pd.DataFrame,
    labels: dict | None = None,
    per_sample: int = 20000,
    ncols: int = 3,
    bins: int = 150,
    seed: int = 0,
) -> plt.Figure:
    """Faceted 2D density with quadrant crosshairs and per-quadrant % labels."""
    labels = _resolve_labels(labels)
    rng = np.random.default_rng(seed)
    n = len(populations)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
    flat = axes.flatten()
    # corner placement (x, y, halign, valign) in axes fraction for each quadrant
    corners = {
        "dn": (0.03, 0.03, "left", "bottom"),
        "x_pos": (0.97, 0.03, "right", "bottom"),
        "dp": (0.97, 0.97, "right", "top"),
        "y_pos": (0.03, 0.97, "left", "top"),
    }
    for ax, p in zip(flat, populations):
        x = xform.apply(_subsample(p["events"][x_channel].to_numpy(dtype=float), per_sample, rng))
        y = xform.apply(_subsample(p["events"][y_channel].to_numpy(dtype=float), per_sample, rng))
        density_plot(
            x, y, ax=ax, xlabel=f"{x_channel} (logicle)", ylabel=f"{y_channel} (logicle)",
            title=p["sample_id"], bins=bins, clip_percentile=None,
        )
        ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
        ax.set_ylim(_LOGICLE_LO, _LOGICLE_HI)
        ax.axvline(x_threshold, color="k", ls="--", lw=0.8)
        ax.axhline(y_threshold, color="k", ls="--", lw=0.8)
        row = quad_df.loc[quad_df["sample_id"] == p["sample_id"]].iloc[0]
        for k in QUADRANT_KEYS:
            fx, fy, ha, va = corners[k]
            ax.text(
                fx, fy, f"{labels[k]}\n{row[f'pct_{labels[k]}']:.1f}%",
                transform=ax.transAxes, ha=ha, va=va, fontsize=7,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.7),
            )
    for ax in flat[n:]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def plot_quadrant_dose_response(
    ax: plt.Axes,
    quad_df: pd.DataFrame,
    dose_col: str,
    labels: dict | None = None,
    logx: bool = True,
) -> plt.Axes:
    """Stacked area of the four quadrant fractions vs dose (finite doses only)."""
    labels = _resolve_labels(labels)
    sub = quad_df[np.isfinite(pd.to_numeric(quad_df[dose_col], errors="coerce"))].copy()
    sub = sub.sort_values(dose_col)
    ys = [sub[f"pct_{labels[k]}"].to_numpy() for k in QUADRANT_KEYS]
    ax.stackplot(sub[dose_col].to_numpy(), *ys, labels=[labels[k] for k in QUADRANT_KEYS])
    if logx:
        ax.set_xscale("symlog")
    ax.set_xlabel(dose_col)
    ax.set_ylabel("% of singlets")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    return ax
