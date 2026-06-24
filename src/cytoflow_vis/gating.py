"""Interactive polygon gating and gate application.

The polygon is drawn once (on a pooled subsample) and applied identically to
every sample. Vertices are captured with matplotlib's ``PolygonSelector`` and:

  * saved to JSON so the same gate can be re-applied without re-clicking, and
  * turned into a FlowKit ``PolygonGate`` so it can be exported to GatingML /
    FlowJo for a fully reproducible record.

Gate membership itself is evaluated with ``matplotlib.path.Path`` for speed.
"""

from __future__ import annotations

import json
from pathlib import Path

import flowkit as fk
import numpy as np
from matplotlib.path import Path as MplPath


def draw_polygon_gate(
    x,
    y,
    x_channel: str = "FSC-A",
    y_channel: str = "SSC-A",
    bins: int = 300,
    initial=None,
) -> list[tuple[float, float]]:
    """Open an interactive window to draw a polygon gate; return its vertices.

    Controls: left-click to add vertices, drag a vertex to move it, ``Esc`` to
    start over, ``Enter`` (or close the window) when finished.

    Requires an interactive matplotlib backend (e.g. QtAgg). Blocks until the
    window is closed.
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import PolygonSelector

    from cytoflow_vis.plotting import density_plot

    fig, ax = plt.subplots(figsize=(8, 8))
    density_plot(
        x,
        y,
        ax=ax,
        xlabel=x_channel,
        ylabel=y_channel,
        title="Draw cell gate: click vertices • Esc resets • Enter when done",
    )

    state: dict = {"verts": list(initial) if initial else []}

    def on_select(verts):
        state["verts"] = [tuple(map(float, v)) for v in verts]

    selector = PolygonSelector(ax, on_select, useblit=True, props=dict(color="red", linewidth=2))
    if initial:
        selector.verts = [tuple(v) for v in initial]

    def on_key(event):
        if event.key == "enter":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show(block=True)

    if len(state["verts"]) < 3:
        raise RuntimeError("A polygon gate needs at least 3 vertices; none captured.")
    return state["verts"]


def save_gate(
    vertices,
    path: str | Path,
    x_channel: str,
    y_channel: str,
    name: str = "cells",
) -> None:
    """Persist a gate's vertices and axes to JSON."""
    payload = {
        "name": name,
        "x_channel": x_channel,
        "y_channel": y_channel,
        "vertices": [list(map(float, v)) for v in vertices],
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def load_gate(path: str | Path) -> dict:
    """Load a gate saved by :func:`save_gate`."""
    return json.loads(Path(path).read_text())


def points_in_polygon(df, x_channel: str, y_channel: str, vertices) -> np.ndarray:
    """Boolean mask: which rows of ``df`` fall inside the polygon."""
    path = MplPath(np.asarray(vertices, dtype=float))
    pts = df[[x_channel, y_channel]].to_numpy()
    return path.contains_points(pts)


def build_flowkit_polygon_gate(
    name: str,
    x_channel: str,
    y_channel: str,
    vertices,
) -> fk.gates.PolygonGate:
    """Construct a FlowKit PolygonGate (for GatingML / FlowJo export)."""
    dims = [fk.Dimension(x_channel), fk.Dimension(y_channel)]
    verts = [list(map(float, v)) for v in vertices]
    return fk.gates.PolygonGate(name, dimensions=dims, vertices=verts)


def seed_populations(samples) -> list[dict]:
    """Turn loaded samples into starting populations (all events, ungated).

    A *population* is a dict with ``sample_id``, ``filename``, ``conditions``
    and an ``events`` DataFrame. Gates consume a population list and return a
    new one with fewer events, so stages chain: cells -> singlets -> ...
    """
    return [
        {
            "sample_id": s.sample_id,
            "filename": s.filename,
            "conditions": s.conditions,
            "events": s.events(),  # all events, not the subsample
        }
        for s in samples
    ]


def apply_gate(populations: list[dict], vertices, x_channel: str, y_channel: str) -> list[dict]:
    """Filter each population to the events inside the polygon.

    Returns a new population list (same identity/conditions, fewer events). Use
    :func:`seed_populations` to create the first list from loaded samples.
    """
    gated = []
    for p in populations:
        df = p["events"]
        mask = points_in_polygon(df, x_channel, y_channel, vertices)
        gated.append({**p, "events": df[mask]})
    return gated


def apply_saved_gates(populations: list[dict], gate_paths) -> list[dict]:
    """Replay a sequence of saved gate JSON files onto the populations.

    Lets downstream steps (e.g. fluorescence analysis) reconstruct a gated
    population from the sample sheet plus the gates drawn earlier, without
    re-clicking. Gates are applied in the order given.
    """
    for gpath in gate_paths:
        gate = load_gate(gpath)
        populations = apply_gate(
            populations, gate["vertices"], gate["x_channel"], gate["y_channel"]
        )
    return populations
