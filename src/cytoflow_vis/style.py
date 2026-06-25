"""Shared plotting style — adapted from the author's ezplot signature look.

A bold, clean, publication aesthetic: dark ink, no top/right spines, thick
spines and ticks, large bold labels, 300-dpi output. The rc values are scaled
for the single-panel (~9 inch) square figures cytoflow-vis produces; pass a
larger ``scale`` for bigger multi-panel figures.
"""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

INK = "#212121"

# The ezplot signature palette (purple-forward), reused for condition colours.
PALETTE = ["#8069EC", "#EA822C", "#2ECC71", "#D143A4", "#F1C40F", "#34495E", "#648FFF"]

# On-brand density ramp: sparse bins fade toward white (so single-event noise
# stays quiet), the body rises through the palette purple, and the dense core
# deepens to a near-ink indigo — keeping good dynamic range on a white ground.
DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "cyto_purple",
    ["#f4f1fc", "#8069EC", "#4b2fa3", "#1a1335"],
)


def rc(scale: float = 1.0) -> dict:
    """Matplotlib rc dict for the signature style, scaled to figure size.

    ``scale=1.0`` is tuned for a single-panel (~9 inch) figure and is bold
    enough to survive projection and downscaling to a journal column; pass
    ``scale<1`` for the smaller panels of a multi-panel facet.
    """
    return {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "text.color": INK,
        "font.size": 23 * scale,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 34 * scale,
        "axes.labelpad": 10 * scale,
        "axes.labelcolor": INK,
        "axes.labelweight": "bold",
        "axes.titlesize": 36 * scale,
        "axes.titleweight": "bold",
        "axes.titlepad": 16 * scale,
        "axes.linewidth": 4 * scale,
        "axes.edgecolor": INK,
        "xtick.bottom": True,
        "ytick.left": True,
        "xtick.major.pad": 10 * scale,
        "ytick.major.pad": 10 * scale,
        "xtick.labelsize": 23 * scale,
        "ytick.labelsize": 23 * scale,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.size": 12 * scale,
        "ytick.major.size": 12 * scale,
        "xtick.major.width": 4 * scale,
        "ytick.major.width": 4 * scale,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": 23 * scale,
        "legend.frameon": False,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
