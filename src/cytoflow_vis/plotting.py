"""Density plots for scatter channels (e.g. SSC-A vs FSC-A).

At 10K-100K events a raw scatter is hopelessly overplotted, so we bin the
events into a 2D histogram and colour by log density. Styling follows the
shared signature look in :mod:`cytoflow_vis.style`.
"""

from __future__ import annotations

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator

from cytoflow_vis.style import INK


def _si_tick(value, _pos=None) -> str:
    """Format a scatter-axis tick the cytometry way: 50K, 150K, 1M, ..."""
    if value == 0:
        return "0"
    av = abs(value)
    if av >= 1e6:
        return f"{value / 1e6:g}M"
    if av >= 1e3:
        return f"{value / 1e3:g}K"
    return f"{value:g}"


def density_plot(
    x,
    y,
    ax: plt.Axes | None = None,
    xlabel: str = "FSC-A",
    ylabel: str = "SSC-A",
    bins: int = 150,
    title: str | None = None,
    cmap="turbo",
    clip_percentile: float | None = 99.9,
    colorbar: bool = False,
) -> plt.Axes:
    """Plot a 2D log-density histogram of ``x`` vs ``y``.

    ``clip_percentile`` trims the axis limits to that percentile so a handful of
    extreme events don't compress the interesting region; set to ``None`` to
    show the full range. Set ``colorbar=True`` to annotate the density scale.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 9))

    if clip_percentile is not None:
        x_hi = np.percentile(x, clip_percentile)
        y_hi = np.percentile(y, clip_percentile)
        x_lo = min(0.0, float(x.min()))
        y_lo = min(0.0, float(y.min()))
        hist_range = [[x_lo, x_hi], [y_lo, y_hi]]
    else:
        hist_range = None

    _, _, _, mesh = ax.hist2d(
        x, y, bins=bins, range=hist_range, cmap=cmap, cmin=1, norm=LogNorm()
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if hist_range is not None:
        ax.set_xlim(hist_range[0])
        ax.set_ylim(hist_range[1])

    # Scatter channels span 0..~262k; use the cytometry-standard K/M tick
    # labels (50K, 150K, ...) rather than a scientific ×10^5 offset.
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(_si_tick))
    ax.yaxis.set_major_formatter(FuncFormatter(_si_tick))
    ax.set_axisbelow(False)

    if colorbar:
        cbar = ax.figure.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label("events / bin", color=INK, weight="bold")
        cbar.outline.set_edgecolor(INK)
        # LogNorm only labels decade majors; over a narrow count range that
        # leaves 2 ticks. Place readable integer increments across the range.
        vmin = max(1.0, float(mesh.norm.vmin or 1.0))
        vmax = float(mesh.norm.vmax or vmin)
        candidates = np.array([1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300,
                               500, 1000, 2000, 5000, 10000], dtype=float)
        ticks = candidates[(candidates >= vmin) & (candidates <= vmax)]
        if ticks.size >= 2:
            cbar.set_ticks(ticks)
            cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)}"))
            cbar.minorticks_off()

    return ax


def overlay_polygon(
    ax: plt.Axes,
    vertices,
    color: str = INK,
    linewidth: float = 4.0,
    markers: bool = True,
    expand_limits: bool = True,
    **kwargs,
) -> plt.Axes:
    """Draw a closed polygon (the gate) on top of an existing plot.

    The line carries a white outline (path effect) so it stays legible over any
    region of the density colourmap; vertices are marked by default. With
    ``expand_limits`` the axis limits grow to keep the whole gate in view (the
    density plot clips to a percentile, so a gate drawn wider than the cluster
    would otherwise be cropped at the frame).
    """
    if vertices is None or len(vertices) < 2:
        return ax
    verts = np.asarray(vertices, dtype=float)
    closed = np.vstack([verts, verts[0]])  # close the loop

    if expand_limits:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        mx = 0.03 * (x1 - x0)
        my = 0.03 * (y1 - y0)
        ax.set_xlim(min(x0, verts[:, 0].min() - mx), max(x1, verts[:, 0].max() + mx))
        ax.set_ylim(min(y0, verts[:, 1].min() - my), max(y1, verts[:, 1].max() + my))
    stroke = [pe.Stroke(linewidth=linewidth + 2.5, foreground="white"), pe.Normal()]
    ax.plot(
        closed[:, 0],
        closed[:, 1],
        color=color,
        linewidth=linewidth,
        solid_joinstyle="round",
        solid_capstyle="round",
        path_effects=stroke,
        zorder=5,
        **kwargs,
    )
    if markers:
        ax.plot(
            verts[:, 0],
            verts[:, 1],
            linestyle="none",
            marker="o",
            markersize=9,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=2,
            zorder=6,
        )
    return ax
