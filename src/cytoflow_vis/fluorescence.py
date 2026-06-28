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
from matplotlib.colors import LogNorm
from matplotlib.ticker import (
    FixedFormatter,
    FixedLocator,
    FuncFormatter,
    SymmetricalLogLocator,
)

from cytoflow_vis.plotting import _si_tick
from cytoflow_vis.style import CATEGORICAL_PALETTE, INK, sequential_colors

DEFAULT_FLUOR_CHANNELS = ["BL1-A", "RL1-A"]

# Logicle display range (a touch below 0 to show the negative population).
_LOGICLE_LO, _LOGICLE_HI = -0.15, 1.05

# Decades labelled on a logicle (biexponential) axis, the cytometry convention.
_LOGICLE_DECADES = [100, 1000, 10000, 100000]

# Facet tick labels and tick marks start from rc(scale<1), so they shrink with
# the panels; enlarge them (relative to the panel's sizes) for legibility when
# printed small or projected. The decades-only logicle axis (no 0 tick) leaves
# the remaining 10^2..10^5 labels evenly spaced, so they take a generous bump.
_FACET_TICK_SCALE = 1.45  # tick *label* size multiplier
_FACET_TICKMARK_LEN_SCALE = 1.7  # tick *mark* length multiplier
_FACET_TICKMARK_WIDTH_SCALE = 1.0  # tick *mark* width multiplier (kept thin so the
# many minor ticks don't crowd)


def _style_facet_ticks(ax) -> None:
    """Enlarge a facet's tick labels and major/minor tick marks together."""
    base_len = plt.rcParams["xtick.major.size"]
    base_w = plt.rcParams["xtick.major.width"]
    ax.tick_params(which="both", labelsize=plt.rcParams["xtick.labelsize"] * _FACET_TICK_SCALE)
    ax.tick_params(which="major", length=base_len * _FACET_TICKMARK_LEN_SCALE,
                   width=base_w * _FACET_TICKMARK_WIDTH_SCALE)
    ax.tick_params(which="minor", length=base_len * _FACET_TICKMARK_LEN_SCALE * 0.55,
                   width=base_w * _FACET_TICKMARK_WIDTH_SCALE * 0.7)


def logicle_axis(ax, xform, which: str = "x", lo: float = _LOGICLE_LO,
                 hi: float = _LOGICLE_HI, label_zero: bool = True) -> None:
    """Relabel a logicle display axis with biexponential decade ticks.

    The events are plotted in logicle *display* units (~[0, 1]); this maps the
    ticks back to raw values (0, 10^2 ... 10^5) at their display positions, with
    minor ticks at the 2-9 subdivisions, so the axis reads like FlowJo/Attune.
    The logicle linear region compresses 0 close to 10^2; on small facets set
    ``label_zero=False`` for a decades-only axis (10^2..10^5) that drops the 0
    tick entirely — a labelless tick there would misread as 10^1.
    """
    decade_raw = list(_LOGICLE_DECADES)
    decade_labels = [f"$10^{int(round(np.log10(d)))}$" for d in decade_raw]
    major_raw = np.array(([0, *decade_raw] if label_zero else decade_raw), dtype=float)
    labels = ["0", *decade_labels] if label_zero else decade_labels
    minor_raw = np.array(
        [k * d for d in [10, *_LOGICLE_DECADES] for k in range(2, 10)], dtype=float
    )
    major_pos = xform.apply(major_raw)
    minor_pos = xform.apply(minor_raw)

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


# --- Replicate handling ------------------------------------------------------
# Samples sharing every condition except the replicate column are replicates of
# one condition. Display plots show a single *representative* replicate (the one
# whose signal is closest to the group median); quantitative plots aggregate
# across replicates. ``DEFAULT_REPLICATE_COL`` is the sample-sheet column.
DEFAULT_REPLICATE_COL = "replicate"


def _varying_conditions(populations: list[dict], replicate_col: str) -> list[str]:
    keys = [k for k in populations[0]["conditions"] if k != replicate_col]
    return [k for k in keys if len({p["conditions"][k] for p in populations}) > 1]


def representative_population(pops: list[dict], channels: list[str], xform) -> dict:
    """The replicate whose median signal is closest to the group median.

    Signal is summarised as the mean over ``channels`` of each replicate's
    median logicle value, so the chosen panel is a typical, not extreme, one.
    """
    if len(pops) == 1:
        return pops[0]
    scores = []
    for p in pops:
        meds = [
            float(np.median(xform.apply(p["events"][ch].to_numpy(dtype=float))))
            for ch in channels
            if ch in p["events"].columns
        ]
        scores.append(float(np.mean(meds)) if meds else 0.0)
    scores = np.asarray(scores)
    return pops[int(np.argmin(np.abs(scores - float(np.median(scores)))))]


def _fmt_value(v) -> str:
    """Render a condition value cleanly: a whole-number float drops its ``.0``."""
    if isinstance(v, float) and np.isfinite(v) and v.is_integer():
        return str(int(v))
    return str(v)


def _auto_ncols(n: int) -> int:
    """Near-square column count for ``n`` facets (4 -> 2x2, 6 -> 3x2)."""
    if n <= 3:
        return n  # a single row reads best for a handful of panels
    return int(np.ceil(np.sqrt(n)))


def _shared_count_norm(hists, vmin=1.0):
    """Global LogNorm over the facet histograms, so a colour means the same
    absolute ``events / bin`` in every panel (panels are subsampled to a common
    N upstream, which keeps the counts comparable)."""
    vmax = max(vmin, max(float(H.max()) for H in hists))
    return LogNorm(vmin=vmin, vmax=vmax), vmax


def _hist_mesh(ax, H, xe, ye, norm, cmap):
    """Draw one raw 2D-histogram density panel (empty bins left white)."""
    return ax.pcolormesh(xe, ye, np.ma.masked_less(H.T, 1), norm=norm, cmap=cmap,
                         rasterized=True)


def _counts_key(fig, axes, mesh, vmax, cmap) -> None:
    """Shared ``events / bin`` colourbar with readable integer ticks."""
    cbar = fig.colorbar(mesh, ax=axes, fraction=0.046, pad=0.02)
    cbar.set_label("events / bin", color=INK, weight="bold")
    cbar.ax.tick_params(labelsize=plt.rcParams["axes.labelsize"] * 0.8)
    cbar.outline.set_edgecolor(INK)
    # LogNorm only labels decade majors; place readable integer increments across
    # the range so the scale is legible even when it spans < 2 decades.
    candidates = np.array([1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500,
                           1000, 2000, 5000, 10000], dtype=float)
    ticks = candidates[(candidates >= 1) & (candidates <= vmax)]
    if ticks.size >= 2:
        cbar.set_ticks(ticks)
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)}"))
        cbar.minorticks_off()


def condition_groups(
    populations: list[dict],
    replicate_col: str = DEFAULT_REPLICATE_COL,
    unit: str | None = None,
) -> list[tuple[str, list[dict]]]:
    """Group populations by condition (all conditions except the replicate col).

    Returns an ordered list of ``(label, replicate_pops)`` — one entry per
    distinct condition, sorted by condition value. With a single varying
    condition the label is its bare value plus an optional ``unit`` (e.g.
    ``"5 mM"``); with several it joins the values.
    """
    varying = _varying_conditions(populations, replicate_col)
    groups: dict = {}
    for p in populations:
        key = tuple(p["conditions"][k] for k in varying)
        groups.setdefault(key, []).append(p)
    out = []
    for key, pops in sorted(groups.items(), key=lambda kv: (kv[0] is None, kv[0])):
        if not varying:
            label = pops[0]["sample_id"]
        elif len(varying) == 1:
            label = f"{_fmt_value(key[0])} {unit}" if unit else _fmt_value(key[0])
        else:
            label = ", ".join(_fmt_value(v) for v in key)
        out.append((label, pops))
    return out


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
    channel_label: str | None = None,
    colors: list | None = None,
    bins: int = 120,
    per_sample: int = 20000,
    threshold: float | None = None,
    overlap: float = 1.7,
    smooth: float = 1.5,
    replicate_mode: str = "representative",
    seed: int = 0,
) -> plt.Axes:
    """Ridgeline of ``channel`` over a biexponential logicle x-axis.

    One offset, filled ridge per ``group_col`` value (or per sample when no
    group is given), coloured from the signature palettes and outlined in ink.
    ``replicate_mode`` decides how replicates of a group combine into its ridge:
    ``"representative"`` (default) uses the representative replicate, ``"pool"``
    concatenates them. Each ridge's tick is its bare condition value (number or
    string); ``group_label`` names the y-axis (e.g. ``"Dose (mM)"``, defaulting
    to ``group_col``). When ``threshold`` is given, the positive gate is drawn
    and each ridge is annotated with its % positive. Expects the ``rc()`` context.
    """
    rng = np.random.default_rng(seed)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 8))

    # Collect populations per ridge (samples sharing a group value = replicates).
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
        reps = grouped[key]
        if replicate_mode == "pool":
            raw = np.concatenate([g["events"][channel].to_numpy(dtype=float) for g in reps])
        else:  # representative replicate for this group
            raw = representative_population(reps, [channel], xform)["events"][channel].to_numpy(dtype=float)
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
            # Match each label to its ridge colour, with a white halo so it still
            # reads over a coloured peak or the white margin.
            ax.text(_LOGICLE_HI - 0.01, baseline + 0.12 * overlap, f"{pct:.0f}%",
                    ha="right", va="bottom", color=color, fontweight="bold",
                    fontsize=plt.rcParams["xtick.labelsize"], zorder=n + 5,
                    path_effects=[pe.Stroke(linewidth=3.5, foreground="white"), pe.Normal()])

    if threshold is not None:
        ax.axvline(threshold, color=INK, ls=(0, (4, 3)), lw=2.5, zorder=n + 4)

    logicle_axis(ax, xform, which="x")
    # Make the logicle minor ticks visible (matplotlib's default is tiny) so the
    # log subdivisions show — they default to a fraction of the major tick.
    ax.tick_params(axis="x", which="minor", color=INK,
                   length=plt.rcParams["xtick.major.size"] * 0.5,
                   width=plt.rcParams["xtick.major.width"] * 0.55)
    ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
    ax.set_ylim(-0.25, (n - 1) * step + overlap + 0.35)
    ax.set_yticks(baselines)
    ax.set_yticklabels([_fmt_value(k) for k in order])
    ax.tick_params(axis="y", which="both", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(channel_label or channel)
    ax.set_ylabel(group_label or group_col or "Sample")
    return ax


def plot_2d_density(
    populations: list[dict],
    x_channel: str,
    y_channel: str,
    xform,
    per_sample: int = 20000,
    ncols: int | None = None,
    bins: int = 150,
    panel: float = 4.5,
    mode: str = "representative",
    replicate_col: str = DEFAULT_REPLICATE_COL,
    unit: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    cmap: str = "viridis",
    seed: int = 0,
) -> plt.Figure:
    """Faceted logicle 2D density.

    ``mode="representative"`` (default) shows one panel per condition, picking
    the representative replicate; ``mode="all"`` is a QC view with one panel per
    sample. Panel titles are the bare condition value plus an optional ``unit``
    (e.g. ``"5 mM"``). Panels are raw 2D-histogram density (the cytometry
    standard) sharing one absolute ``events / bin`` colour scale, biexponential
    logicle axes, and shared figure-level channel labels; only the outer panels
    carry tick labels. ``ncols`` defaults to a near-square tiling. Expects the
    ``rc()`` context (reduced ``scale``).
    """
    rng = np.random.default_rng(seed)
    if mode == "all":
        panels = [(p["sample_id"], p) for p in populations]
    else:
        panels = [
            (label, representative_population(pops, [x_channel, y_channel], xform))
            for label, pops in condition_groups(populations, replicate_col, unit=unit)
        ]
    n = len(panels)
    ncols = min(ncols or _auto_ncols(n), n)
    nrows = int(np.ceil(n / ncols))
    hist_range = [[_LOGICLE_LO, _LOGICLE_HI], [_LOGICLE_LO, _LOGICLE_HI]]
    binned = []
    for label, p in panels:
        x = xform.apply(_subsample(p["events"][x_channel].to_numpy(dtype=float), per_sample, rng))
        y = xform.apply(_subsample(p["events"][y_channel].to_numpy(dtype=float), per_sample, rng))
        H, xe, ye = np.histogram2d(x, y, bins=bins, range=hist_range)
        binned.append((label, H, xe, ye))
    norm, vmax = _shared_count_norm([H for _, H, _, _ in binned])

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(panel * ncols, panel * nrows),
        squeeze=False, sharex=True, sharey=True, layout="constrained",
    )
    fig.get_layout_engine().set(hspace=0.10)  # breathing room between rows
    flat = axes.flatten()
    mesh = None
    for ax, (label, H, xe, ye) in zip(flat, binned):
        mesh = _hist_mesh(ax, H, xe, ye, norm, cmap)
        ax.set_title(label, pad=3)  # sits right atop its panel
        ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
        ax.set_ylim(_LOGICLE_LO, _LOGICLE_HI)
        logicle_axis(ax, xform, which="x", label_zero=False)
        logicle_axis(ax, xform, which="y", label_zero=False)
        _style_facet_ticks(ax)
        ax.label_outer()  # only outer panels keep tick labels
    for ax in flat[n:]:
        ax.axis("off")
    label_kw = dict(fontsize=plt.rcParams["axes.labelsize"], fontweight="bold")
    fig.supxlabel(xlabel or x_channel, **label_kw)
    fig.supylabel(ylabel or y_channel, **label_kw)
    _counts_key(fig, axes, mesh, vmax, cmap)
    return fig


def _dose_linthresh(values) -> float:
    """Symlog linear-threshold for a dose axis = the smallest positive dose.

    Putting the linear region just below the lowest real dose means every
    positive dose sits in the log region (properly spread), with 0 alone in the
    narrow linear zone — instead of fractional doses bunching against 0.
    """
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy()
    pos = v[np.isfinite(v) & (v > 0)]
    return float(pos.min()) if pos.size else 1.0


def _symlog_dose_axis(ax, linthresh: float | None = None) -> None:
    """Put the dose x-axis on a symlog scale with 2-9 minor ticks per decade.

    Symlog (not log) so a zero dose still has a place; ``linthresh`` sets where
    the log region begins (default the smallest positive dose), so fractional
    doses spread out rather than bunch against 0. The minor ticks mark the log
    subdivisions, signalling the scale and helping read intermediate doses.
    """
    ax.set_xscale("symlog", linthresh=linthresh) if linthresh else ax.set_xscale("symlog")
    ax.xaxis.set_minor_locator(
        SymmetricalLogLocator(ax.xaxis.get_transform(), subs=np.arange(2, 10))
    )
    major = plt.rcParams["xtick.major.size"]
    ax.tick_params(axis="x", which="minor",
                   length=major * 0.68, width=plt.rcParams["xtick.major.width"] * 0.72,
                   color=INK)
    # Hold the left edge just inside -linthresh: 0 keeps a little margin, but the
    # symmetric negative decade tick (e.g. -10^-1) — meaningless for a dose — is
    # pushed out of view.
    lt = linthresh or 1.0
    ax.set_xlim(left=-0.5 * lt)


def _dose_summary(sub: pd.DataFrame, dose_col: str, y_col: str) -> pd.DataFrame:
    """Per-dose mean, SD and replicate count of ``y_col`` (sorted by dose)."""
    x = pd.to_numeric(sub[dose_col], errors="coerce")
    g = sub.assign(_dose=x).groupby("_dose")[y_col]
    agg = g.agg(mean="mean", sd="std", n="count").reset_index().sort_values("_dose")
    agg["sd"] = agg["sd"].fillna(0.0)  # a single replicate has no spread
    return agg


def _jitter_dose(x: np.ndarray, linthresh: float, spread: float = 0.05) -> np.ndarray:
    """Spread co-located replicate x-values apart so none hides behind another.

    Replicates at the same dose are dodged symmetrically: multiplicatively for
    non-zero doses (even spacing on the log axis) and additively at zero (which
    sits in the symlog linear region). The zero dodge is scaled to ``linthresh``
    so it stays a small fraction of the linear region — never reaching the
    smallest positive dose. A lone point is left untouched.
    """
    x = np.asarray(x, dtype=float)
    out = x.copy()
    for v in np.unique(x):
        idx = np.where(x == v)[0]
        if idx.size <= 1:
            continue
        offs = np.linspace(-1.0, 1.0, idx.size)
        out[idx] = v * 10.0 ** (spread * offs) if v != 0 else 0.11 * linthresh * offs
    return out


# Distinct marker shapes per series, so overlapping lines stay distinguishable
# by shape (redundant with colour) — robust for colour-blind and grayscale view.
_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


def _dose_series(ax, df, dose_col, y_col, color, linthresh, marker="o", label=None) -> set:
    """Draw one dose-response series: faint replicate points + mean ± SD curve.

    Replicate points are jittered horizontally so a tight cluster stays visible
    beside (not behind) the mean marker; ``marker`` shapes the series so it reads
    apart from others even where lines overlap. Returns the set of replicate
    counts seen (for an ``n`` annotation).
    """
    x_pts = pd.to_numeric(df[dose_col], errors="coerce").to_numpy()
    ax.scatter(_jitter_dose(x_pts, linthresh), df[y_col], s=42, color=color, alpha=0.55,
               marker=marker, linewidths=0, zorder=3)
    agg = _dose_summary(df, dose_col, y_col)
    # A white underline gives the line its halo (so crossing lines stay separable)
    # without haloing the markers. The coloured line *with* its marker then carries
    # the label, so the legend handle is a proper line+marker (not marker-only).
    ax.plot(agg["_dose"], agg["mean"], color="white", lw=6.5, zorder=3.8,
            solid_capstyle="round", solid_joinstyle="round")
    ax.errorbar(agg["_dose"], agg["mean"], yerr=agg["sd"], fmt="none", ecolor=color,
                elinewidth=2.5, capsize=6, capthick=2.5, zorder=3.9)
    ax.plot(agg["_dose"], agg["mean"], color=color, lw=3.5, marker=marker, markersize=11,
            markeredgecolor="white", markeredgewidth=1.0, zorder=4, label=label)
    return {int(v) for v in agg["n"]}


def plot_dose_response(
    ax: plt.Axes,
    stats_df: pd.DataFrame,
    dose_col: str,
    y_col: str,
    group_col: str | None = None,
    group_label: str | None = None,
    colors: list | None = None,
    dose_label: str | None = None,
    y_label: str | None = None,
    percent: bool = False,
    show_n: bool = True,
    show_legend: bool = True,
    logx: bool = True,
) -> plt.Axes:
    """Dose-response of ``y_col`` vs ``dose_col`` with per-replicate detail.

    Each biological replicate is drawn as a faint individual point; the mean
    across replicates is the bold connected curve with mean ± SD error bars (no
    pseudoreplication — replicates are shown, not pooled). One series per
    ``group_col`` value (e.g. cell line); a lone series uses a single accent.
    ``group_label`` titles the legend (defaulting to ``group_col``). Colours
    follow the data-type convention (categorical groups -> Okabe-Ito);
    ``percent`` pins the y-axis to 0-100 %, otherwise raw MFI gets K/M ticks.
    Expects the ``rc()`` context.
    """
    if group_col and group_col in stats_df and stats_df[group_col].nunique() > 1:
        series = [(gval, sub) for gval, sub in stats_df.groupby(group_col)]
    else:
        series = [(None, stats_df)]

    if colors is None:
        colors = (
            [CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)] for i in range(len(series))]
            if len(series) > 1
            else [CATEGORICAL_PALETTE[0]]  # Okabe deep blue accent
        )

    linthresh = _dose_linthresh(stats_df[dose_col])
    multi = len(series) > 1
    n_reps = set()
    for i, (gval, sub) in enumerate(series):
        n_reps |= _dose_series(ax, sub, dose_col, y_col, colors[i], linthresh,
                               marker=_MARKERS[i % len(_MARKERS)] if multi else "o",
                               label=str(gval) if gval is not None else None)

    if logx:
        # log region starts at the smallest positive dose so 0/fractional doses
        # don't bunch up; symlog still gives a zero dose a place.
        _symlog_dose_axis(ax, linthresh)
    if percent:
        ax.set_ylim(-4, 104)  # headroom so markers/line at 0 % and 100 % aren't clipped
    else:
        ax.set_ylim(bottom=0)  # anchor at 0 so a flat channel reads as flat
        ax.yaxis.set_major_formatter(FuncFormatter(_si_tick))
    ax.set_xlabel(dose_label or dose_col)
    ax.set_ylabel(y_label or y_col)

    # n is the replicate count; surface it when uniform across doses/groups. In a
    # grid the caller reports it once at the figure level instead (show_n=False),
    # and likewise draws one shared group legend (show_legend=False).
    n_note = f"n = {n_reps.pop()}" if len(n_reps) == 1 else None
    if show_legend and len(series) > 1:
        gl = group_label or group_col
        title = f"{gl}  ({n_note})" if show_n and n_note else gl
        ax.legend(title=title, loc="best").get_title().set_fontweight("bold")
    elif show_legend and show_n and n_note:
        # A frameless single-entry legend so "best" placement dodges the curve.
        (handle,) = ax.plot([], [], linestyle="none")
        ax.legend([handle], [n_note], loc="best", handlelength=0, handletextpad=0,
                  borderaxespad=0.8)
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
    ncols: int | None = None,
    bins: int = 150,
    panel: float = 4.5,
    replicate_col: str = DEFAULT_REPLICATE_COL,
    unit: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    cmap: str = "viridis",
    seed: int = 0,
) -> plt.Figure:
    """Faceted 2D density with quadrant crosshairs and per-quadrant % labels.

    One panel per condition (the representative replicate, matching the 2D
    density facets); the crosshairs mark the two thresholds and each corner is
    annotated with that panel's quadrant percentage. Panels are raw 2D-histogram
    density sharing one absolute ``events / bin`` scale, with shared biexponential
    logicle axes, channel-biology figure labels, and a near-square tiling.
    Expects ``rc()``.
    """
    labels = _resolve_labels(labels)
    rng = np.random.default_rng(seed)
    panels = [
        (label, representative_population(pops, [x_channel, y_channel], xform))
        for label, pops in condition_groups(populations, replicate_col, unit=unit)
    ]
    n = len(panels)
    ncols = min(ncols or _auto_ncols(n), n)
    nrows = int(np.ceil(n / ncols))
    hist_range = [[_LOGICLE_LO, _LOGICLE_HI], [_LOGICLE_LO, _LOGICLE_HI]]
    binned = []
    for label, p in panels:
        x = xform.apply(_subsample(p["events"][x_channel].to_numpy(dtype=float), per_sample, rng))
        y = xform.apply(_subsample(p["events"][y_channel].to_numpy(dtype=float), per_sample, rng))
        H, xe, ye = np.histogram2d(x, y, bins=bins, range=hist_range)
        binned.append((label, p, H, xe, ye))
    norm, vmax = _shared_count_norm([H for _, _, H, _, _ in binned])

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(panel * ncols, panel * nrows),
        squeeze=False, sharex=True, sharey=True, layout="constrained",
    )
    fig.get_layout_engine().set(hspace=0.10)
    flat = axes.flatten()
    # corner placement (x, y, halign, valign) in axes fraction for each quadrant
    corners = {
        "dn": (0.04, 0.04, "left", "bottom"),
        "x_pos": (0.96, 0.04, "right", "bottom"),
        "dp": (0.96, 0.96, "right", "top"),
        "y_pos": (0.04, 0.96, "left", "top"),
    }
    mesh = None
    for ax, (label, p, H, xe, ye) in zip(flat, binned):
        mesh = _hist_mesh(ax, H, xe, ye, norm, cmap)
        ax.axvline(x_threshold, color=INK, ls=(0, (4, 3)), lw=2.0, zorder=4)
        ax.axhline(y_threshold, color=INK, ls=(0, (4, 3)), lw=2.0, zorder=4)
        row = quad_df.loc[quad_df["sample_id"] == p["sample_id"]].iloc[0]
        for k in QUADRANT_KEYS:
            fx, fy, ha, va = corners[k]
            ax.text(
                fx, fy, f"{labels[k]}\n{row[f'pct_{labels[k]}']:.0f}%",
                transform=ax.transAxes, ha=ha, va=va, color=INK, fontweight="bold",
                fontsize=plt.rcParams["xtick.labelsize"], zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=INK, lw=1.2, alpha=0.85),
            )
        ax.set_title(label, pad=3)
        ax.set_xlim(_LOGICLE_LO, _LOGICLE_HI)
        ax.set_ylim(_LOGICLE_LO, _LOGICLE_HI)
        logicle_axis(ax, xform, which="x", label_zero=False)
        logicle_axis(ax, xform, which="y", label_zero=False)
        _style_facet_ticks(ax)
        ax.label_outer()
    for ax in flat[n:]:
        ax.axis("off")
    label_kw = dict(fontsize=plt.rcParams["axes.labelsize"], fontweight="bold")
    fig.supxlabel(xlabel or x_channel, **label_kw)
    fig.supylabel(ylabel or y_channel, **label_kw)
    _counts_key(fig, axes, mesh, vmax, cmap)
    return fig


def plot_quadrant_dose_response(
    ax: plt.Axes,
    quad_df: pd.DataFrame,
    dose_col: str,
    labels: dict | None = None,
    colors: list | None = None,
    dose_label: str | None = None,
    legend_title: str = "Population",
    show_n: bool = True,
    logx: bool = True,
) -> plt.Axes:
    """Per-quadrant % vs dose: one curve per quadrant with replicate detail.

    Each of the four quadrant fractions is a series with faint per-replicate
    points and a mean ± SD curve (the same construction as the dose-response),
    coloured from the categorical palette. The legend is titled ``legend_title``;
    with ``show_n`` the replicate count is appended (``"Population (n = 3)"``),
    else the caller reports n once at the figure level. Expects the ``rc()`` context.
    """
    labels = _resolve_labels(labels)
    sub = quad_df[np.isfinite(pd.to_numeric(quad_df[dose_col], errors="coerce"))].copy()
    if colors is None:
        colors = [CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)] for i in range(len(QUADRANT_KEYS))]
    linthresh = _dose_linthresh(sub[dose_col])
    n_reps = set()
    for i, k in enumerate(QUADRANT_KEYS):
        n_reps |= _dose_series(ax, sub, dose_col, f"pct_{labels[k]}", colors[i], linthresh,
                               marker=_MARKERS[i % len(_MARKERS)], label=labels[k])
    if logx:
        _symlog_dose_axis(ax, linthresh)
    ax.set_xlabel(dose_label or dose_col)
    ax.set_ylabel("% of singlets")
    ax.set_ylim(-4, 104)  # headroom so a quadrant at 0 % or 100 % isn't clipped
    title = f"{legend_title} (n = {n_reps.pop()})" if show_n and len(n_reps) == 1 else legend_title
    leg = ax.legend(title=title, loc="center left", bbox_to_anchor=(1.02, 0.5))
    leg.get_title().set_fontweight("bold")
    return ax
