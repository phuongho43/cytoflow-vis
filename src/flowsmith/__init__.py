"""flowsmith: interactive gating and visualization for flow cytometry FCS files."""

from flowsmith.io import LoadedSample, load_sample_sheet, load_samples
from flowsmith import analysis, fluorescence, spillover
from flowsmith.analysis import AnalysisContext, REGISTRY, register, run_analysis
from flowsmith.spillover import compute_spillover_matrix
from flowsmith.gating import (
    apply_gate,
    apply_saved_gates,
    build_flowkit_polygon_gate,
    draw_polygon_gate,
    load_gate,
    points_in_polygon,
    save_gate,
    seed_populations,
)
from flowsmith.plotting import density_plot, overlay_polygon

__all__ = [
    "LoadedSample",
    "load_sample_sheet",
    "load_samples",
    "analysis",
    "fluorescence",
    "spillover",
    "AnalysisContext",
    "REGISTRY",
    "register",
    "run_analysis",
    "compute_spillover_matrix",
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
