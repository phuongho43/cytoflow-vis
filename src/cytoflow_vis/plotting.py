"""Density plots for scatter channels (e.g. SSC-A vs FSC-A).

At 10K-100K events a raw scatter is hopelessly overplotted, so we bin the
events into a 2D histogram and colour by log density.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def density_plot(
    x,
    y,
    ax: plt.Axes | None = None,
    xlabel: str = "FSC-A",
    ylabel: str = "SSC-A",
    bins: int = 300,
    title: str | None = None,
    cmap: str = "turbo",
    clip_percentile: float | None = 99.9,
) -> plt.Axes:
    """Plot a 2D log-density histogram of ``x`` vs ``y``.

    ``clip_percentile`` trims the axis limits to that percentile so a handful of
    extreme events don't compress the interesting region; set to ``None`` to
    show the full range.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    if clip_percentile is not None:
        x_hi = np.percentile(x, clip_percentile)
        y_hi = np.percentile(y, clip_percentile)
        x_lo = min(0.0, float(x.min()))
        y_lo = min(0.0, float(y.min()))
        hist_range = [[x_lo, x_hi], [y_lo, y_hi]]
    else:
        hist_range = None

    ax.hist2d(x, y, bins=bins, range=hist_range, cmap=cmap, cmin=1, norm=LogNorm())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if hist_range is not None:
        ax.set_xlim(hist_range[0])
        ax.set_ylim(hist_range[1])
    return ax


def overlay_polygon(ax: plt.Axes, vertices, color: str = "red", **kwargs) -> plt.Axes:
    """Draw a closed polygon (the gate) on top of an existing plot."""
    if vertices is None or len(vertices) < 2:
        return ax
    verts = np.asarray(vertices, dtype=float)
    closed = np.vstack([verts, verts[0]])  # close the loop
    ax.plot(closed[:, 0], closed[:, 1], color=color, linewidth=2, **kwargs)
    return ax
