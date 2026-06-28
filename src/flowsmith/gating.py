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
from scipy.stats import chi2


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

    from flowsmith.plotting import density_plot

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


# -- automatic gating ---------------------------------------------------------

def _robust_gaussian(points: np.ndarray, support_frac: float = 0.5, iters: int = 10):
    """Robust mean & covariance of the dense core (MCD-style concentration steps).

    Iteratively keeps the ``support_frac`` of points closest (Mahalanobis) to the
    current fit and refits, so a separate debris cloud or scattered outliers do
    not drag the estimate. The covariance is then rescaled by the standard
    consistency factor so it matches a full Gaussian rather than the trimmed core.
    """
    P = np.asarray(points, dtype=float)
    n = len(P)
    h = max(int(support_frac * n), P.shape[1] + 1)
    mu = np.median(P, axis=0)
    mad = np.median(np.abs(P - mu), axis=0) * 1.4826
    mad[mad == 0] = 1.0
    cov = np.diag(mad ** 2)
    for _ in range(iters):
        diff = P - mu
        d2 = np.einsum("ij,jk,ik->i", diff, np.linalg.pinv(cov), diff)
        core = P[np.argsort(d2)[:h]]
        mu = core.mean(axis=0)
        cov = np.cov(core.T)
    # Consistency correction: the trimmed core underestimates spread.
    q = chi2.ppf(h / n, 2)
    cov *= (h / n) / chi2.cdf(q, 4)
    return mu, cov


def auto_ellipse_gate(x, y, coverage: float = 0.95, n_vertices: int = 40,
                      support_frac: float = 0.5) -> list[tuple[float, float]]:
    """Ellipse polygon around the dense main population (e.g. cells vs debris).

    Fits a robust 2D Gaussian to the densest cluster and returns the iso-density
    ellipse enclosing ``coverage`` of it as polygon vertices — a hands-off
    replacement for hand-drawing the FSC/SSC cell gate.
    """
    P = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    mu, cov = _robust_gaussian(P, support_frac=support_frac)
    vals, vecs = np.linalg.eigh(cov)
    vals = np.clip(vals, 1e-12, None)
    r = np.sqrt(chi2.ppf(coverage, 2))
    t = np.linspace(0, 2 * np.pi, n_vertices, endpoint=False)
    circle = np.column_stack([np.cos(t), np.sin(t)])
    ellipse = (circle * (r * np.sqrt(vals))) @ vecs.T + mu
    return [tuple(map(float, v)) for v in ellipse]


def auto_singlet_gate(x, y, k: float = 4.0, x_pct=(0.5, 99.5)) -> list[tuple[float, float]]:
    """Constant-width band around the FSC-H/FSC-A singlet diagonal.

    Singlets fall on a line ``FSC-H ≈ slope·FSC-A`` with roughly constant scatter;
    doublets sit well below it. The slope is the **median height/area ratio**
    (robust to a doublet minority, unlike a least-squares fit they would tilt),
    and events are kept within ``k`` robust SDs of that line — a parallelogram of
    constant vertical width, not an origin wedge that flares out at high FSC-A.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = x > 0
    x, y = x[ok], y[ok]

    slope = float(np.median(y / x))  # robust height/area ratio
    resid = y - slope * x
    c = np.median(resid)  # robust offset of the singlet line
    s = np.median(np.abs(resid - c)) * 1.4826 or 1.0  # robust SD about the line
    on_line = np.abs(resid - c) < 3 * s  # the singlets (drops off-line doublets)
    x_lo, x_hi = np.percentile(x[on_line] if on_line.sum() >= 10 else x, list(x_pct))

    def line(xx):
        return slope * xx + c

    return [(float(x_lo), float(line(x_lo) - k * s)), (float(x_hi), float(line(x_hi) - k * s)),
            (float(x_hi), float(line(x_hi) + k * s)), (float(x_lo), float(line(x_lo) + k * s))]


def auto_gate(stage_name: str, x, y) -> list[tuple[float, float]]:
    """Compute an automatic polygon gate for a standard stage.

    The ``singlets`` stage uses the height/area ratio wedge; every other stage
    (e.g. ``cells``) uses the robust density ellipse.
    """
    if stage_name == "singlets":
        return auto_singlet_gate(x, y)
    return auto_ellipse_gate(x, y)


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
