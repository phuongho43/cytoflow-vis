"""cytoflow-vis: interactive gating and visualization for flow cytometry FCS files."""

from cytoflow_vis.io import LoadedSample, load_sample_sheet, load_samples
from cytoflow_vis import fluorescence
from cytoflow_vis.gating import (
    apply_gate,
    apply_saved_gates,
    build_flowkit_polygon_gate,
    draw_polygon_gate,
    load_gate,
    points_in_polygon,
    save_gate,
    seed_populations,
)
from cytoflow_vis.plotting import density_plot, overlay_polygon

__all__ = [
    "LoadedSample",
    "load_sample_sheet",
    "load_samples",
    "fluorescence",
    "apply_gate",
    "apply_saved_gates",
    "build_flowkit_polygon_gate",
    "draw_polygon_gate",
    "load_gate",
    "points_in_polygon",
    "save_gate",
    "seed_populations",
    "density_plot",
    "overlay_polygon",
]
