"""Composable, config-selectable fluorescence analyses.

Each analysis is a small function registered under a ``kind`` name. A config
file lists the analyses to run (with parameters); the runner builds an
:class:`AnalysisContext` from the gated populations and dispatches each one.

Add a new analysis by writing a function decorated with ``@register("name")``
that takes the context plus keyword parameters and writes its outputs to
``ctx.out_dir``. It then becomes available to every experiment's config.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from cytoflow_vis import fluorescence as fl
from cytoflow_vis.style import rc


@dataclass
class AnalysisContext:
    """Shared inputs handed to every analysis."""

    populations: list[dict]  # gated singlets, with .conditions + .events
    out_dir: Path
    xform: object  # logicle transform
    channels: list[str]  # default fluorescence channels of interest
    control_id: str | None  # default negative control sample id
    group_col: str | None  # default condition to colour/group by
    dose_col: str | None  # default numeric condition for dose-response
    per_sample: int = 20000
    positive_percentile: float = 99.0


# -- registry -----------------------------------------------------------------

REGISTRY: dict[str, callable] = {}


def register(kind: str):
    def deco(fn):
        REGISTRY[kind] = fn
        return fn
    return deco


def run_analysis(kind: str, ctx: AnalysisContext, params: dict) -> str:
    if kind not in REGISTRY:
        raise ValueError(f"unknown analysis kind {kind!r}; available: {sorted(REGISTRY)}")
    return REGISTRY[kind](ctx, **params)


def _safe(name: str) -> str:
    return name.replace("/", "_")


# -- column auto-detection ----------------------------------------------------

def _is_numeric(values) -> bool:
    try:
        [float(v) for v in values]
        return True
    except (TypeError, ValueError):
        return False


def detect_columns(populations, group_arg=None, dose_arg=None):
    """Auto-pick dose (numeric, varying) and group (categorical, varying) columns."""
    cond_keys = list(populations[0]["conditions"].keys())

    def uniques(key):
        return {p["conditions"][key] for p in populations}

    varying = [k for k in cond_keys if len(uniques(k)) > 1]

    dose = dose_arg
    if dose is None:
        numeric_varying = [k for k in varying if _is_numeric(uniques(k))]
        if "concentration" in numeric_varying:
            dose = "concentration"
        elif numeric_varying:
            dose = numeric_varying[0]

    group = group_arg
    if group is None:
        candidates = [k for k in varying if k != dose and not _is_numeric(uniques(k))]
        group = candidates[0] if candidates else None
    return group, dose


def detect_control(populations, control_arg=None, dose_col=None):
    if control_arg is not None:
        return control_arg
    if dose_col is not None:
        return min(populations, key=lambda p: float(p["conditions"][dose_col]))["sample_id"]
    return populations[0]["sample_id"]


def dose_determines_group(populations, group_col, dose_col) -> bool:
    """True if each dose value maps to a single group value (group redundant)."""
    if group_col is None or dose_col is None:
        return False
    mapping = {}
    for p in populations:
        dose, grp = p["conditions"][dose_col], p["conditions"][group_col]
        if dose in mapping and mapping[dose] != grp:
            return False
        mapping[dose] = grp
    return True


# -- analyses -----------------------------------------------------------------

@register("mfi_pct")
def _mfi_pct(ctx, channels=None, control=None, positive_percentile=None, out="fluor_stats.csv"):
    """Per-sample MFI and % positive table."""
    channels = channels or ctx.channels
    control = control if control is not None else ctx.control_id
    pct = positive_percentile if positive_percentile is not None else ctx.positive_percentile
    stats, _ = fl.compute_stats(ctx.populations, channels, ctx.xform, control_id=control, positive_percentile=pct)
    path = ctx.out_dir / out
    stats.to_csv(path, index=False)
    return f"mfi_pct -> {path.name}"


@register("histograms")
def _histograms(ctx, channels=None, group=None, group_label=None, colors=None,
                control=None, positive_percentile=None, per_sample=None):
    """One logicle ridgeline PNG per channel, one ridge per condition value.

    ``group`` is the condition column to split on (auto-detected if omitted);
    ``group_label`` sets the y-axis title with units (e.g. ``"Dose (mM)"``).
    Colours auto-pick by data type (numeric -> sequential ramp, string ->
    categorical palette); ``colors`` overrides with an explicit list.
    """
    channels = channels or ctx.channels
    group = group if group is not None else ctx.group_col
    control = control if control is not None else ctx.control_id
    pct = positive_percentile if positive_percentile is not None else ctx.positive_percentile
    per_sample = per_sample or ctx.per_sample
    thresholds = (
        fl.control_thresholds(ctx.populations, channels, ctx.xform, control, pct) if control else {}
    )
    n = fl.n_groups(ctx.populations, group)
    height = max(4.0, 1.05 * n + 1.8)  # ridgeline grows with the number of rows
    written = []
    for ch in channels:
        with mpl.rc_context(rc()):
            fig, ax = plt.subplots(figsize=(9, height))
            fl.plot_histograms(
                ctx.populations, ch, ctx.xform, ax=ax, group_col=group,
                group_label=group_label, colors=colors, per_sample=per_sample,
                threshold=thresholds.get(ch),
            )
            path = ctx.out_dir / f"hist_{_safe(ch)}.png"
            fig.savefig(path)
            plt.close(fig)
        written.append(path.name)
    return "histograms -> " + ", ".join(written)


@register("density_2d")
def _density_2d(ctx, x=None, y=None, per_sample=None, out="density_2d.png"):
    """Faceted 2D logicle density, one panel per sample."""
    x = x or ctx.channels[0]
    y = y or (ctx.channels[1] if len(ctx.channels) > 1 else ctx.channels[0])
    per_sample = per_sample or ctx.per_sample
    fig = fl.plot_2d_density(ctx.populations, x, y, ctx.xform, per_sample=per_sample)
    path = ctx.out_dir / out
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return f"density_2d -> {path.name}"


@register("dose_response")
def _dose_response(ctx, channels=None, dose=None, group=None, control=None,
                   positive_percentile=None, out="dose_response.png"):
    """MFI and % positive vs dose, one connected curve through the control."""
    channels = channels or ctx.channels
    dose = dose if dose is not None else ctx.dose_col
    if dose is None:
        return "dose_response -> skipped (no dose column)"
    group = group if group is not None else ctx.group_col
    control = control if control is not None else ctx.control_id
    pct = positive_percentile if positive_percentile is not None else ctx.positive_percentile
    stats, thresholds = fl.compute_stats(
        ctx.populations, channels, ctx.xform, control_id=control, positive_percentile=pct
    )
    dr_group = None if dose_determines_group(ctx.populations, group, dose) else group
    metrics = [("MFI", "MFI")] + ([("pct_pos", "% positive")] if thresholds else [])
    fig, axes = plt.subplots(len(metrics), len(channels),
                             figsize=(6 * len(channels), 5 * len(metrics)), squeeze=False)
    for r, (prefix, ylabel) in enumerate(metrics):
        for c, ch in enumerate(channels):
            ax = axes[r][c]
            fl.plot_dose_response(ax, stats, dose, f"{prefix}_{ch}", group_col=dr_group)
            ax.set_ylabel(f"{ylabel} {ch}")
    fig.tight_layout()
    path = ctx.out_dir / out
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return f"dose_response -> {path.name}"


@register("quadrant")
def _quadrant(ctx, x, y, control=None, positive_percentile=None, labels=None,
              dose=None, out="quadrant.csv"):
    """Two-threshold quadrant analysis: % in each of four populations.

    Writes a per-sample CSV, a faceted density plot with crosshairs and
    per-quadrant labels, and (if a dose column is available) a stacked
    dose-response of the quadrant fractions.
    """
    control = control if control is not None else ctx.control_id
    pct = positive_percentile if positive_percentile is not None else ctx.positive_percentile
    thr = fl.control_thresholds(ctx.populations, [x, y], ctx.xform, control, pct)
    quad = fl.quadrant_stats(ctx.populations, x, y, ctx.xform, thr[x], thr[y], labels=labels)
    csv_path = ctx.out_dir / out
    quad.to_csv(csv_path, index=False)
    written = [csv_path.name]

    fig = fl.plot_quadrants(
        ctx.populations, x, y, ctx.xform, thr[x], thr[y], quad, labels=labels,
        per_sample=ctx.per_sample,
    )
    plot_path = ctx.out_dir / f"{csv_path.stem}_density.png"
    fig.savefig(plot_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    written.append(plot_path.name)

    dose = dose if dose is not None else ctx.dose_col
    if dose is not None:
        fig, ax = plt.subplots(figsize=(8, 6))
        fl.plot_quadrant_dose_response(ax, quad, dose, labels=labels)
        dr_path = ctx.out_dir / f"{csv_path.stem}_dose_response.png"
        fig.savefig(dr_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        written.append(dr_path.name)
    return "quadrant -> " + ", ".join(written)
