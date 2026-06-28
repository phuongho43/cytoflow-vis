"""Standalone entry point: draw a sequential gating hierarchy and apply it.

Usage (after `uv sync`):

    uv run gate-cells --sheet samples.csv --data data/ --out results/

Gating stages (each drawn once on a pooled subsample, then applied to every
sample's full event set):

    1. cells    SSC-A vs FSC-A   -- exclude low-scatter debris
    2. singlets FSC-H  vs FSC-A  -- keep the diagonal, exclude doublets/aggregates

Each stage's polygon is filtered from the *previous* stage's population, so the
singlet gate is drawn on the cells only. Drawn gates are saved per stage as
``<stage>_gate.json`` and reused on the next run; pass --redraw to draw fresh.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Stage:
    """One gate in the hierarchy."""

    name: str
    x: str
    y: str


DEFAULT_STAGES = [
    Stage("cells", "FSC-A", "SSC-A"),
    Stage("singlets", "FSC-A", "FSC-H"),
]


def pool_events(populations, x_channel, y_channel, per_sample=5000, seed=0):
    """Concatenate a random subsample of (x, y) events across all populations."""
    rng = np.random.default_rng(seed)
    xs, ys = [], []
    for p in populations:
        df = p["events"]
        if len(df) == 0:
            continue
        n = min(per_sample, len(df))
        idx = rng.choice(len(df), size=n, replace=False)
        xs.append(df[x_channel].to_numpy()[idx])
        ys.append(df[y_channel].to_numpy()[idx])
    return np.concatenate(xs), np.concatenate(ys)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Sequential FSC/SSC gating across many samples.")
    p.add_argument("--sheet", required=True, help="Sample sheet CSV (needs a 'filename' column).")
    p.add_argument("--data", required=True, help="Directory containing the FCS files.")
    p.add_argument("--out", default="results", help="Output directory (default: results/).")
    p.add_argument("--gates", default=None, help="Directory for <stage>_gate.json files (default: --out).")
    p.add_argument("--redraw", action="store_true", help="Re-draw every gate even if its file exists.")
    p.add_argument("--per-sample", type=int, default=5000, help="Events per file pooled for each draw plot.")
    p.add_argument("--subsample", type=int, default=20000, help="Events FlowKit stores per file.")
    p.add_argument("--no-save-events", action="store_true", help="Skip writing the final per-sample CSVs.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    stages = DEFAULT_STAGES

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    gates_dir = Path(args.gates) if args.gates else out_dir
    gates_dir.mkdir(parents=True, exist_ok=True)

    def gate_file(stage: Stage) -> Path:
        return gates_dir / f"{stage.name}_gate.json"

    # If any stage needs drawing, switch to the interactive Qt backend *before*
    # pyplot is imported anywhere (you can't change backend afterwards).
    need_draw = any(args.redraw or not gate_file(s).exists() for s in stages)
    import matplotlib

    if need_draw:
        matplotlib.use("QtAgg")

    from cytoflow_vis.gating import (
        apply_gate,
        build_flowkit_polygon_gate,
        draw_polygon_gate,
        load_gate,
        save_gate,
        seed_populations,
    )
    from cytoflow_vis.io import load_samples
    from cytoflow_vis.plotting import density_plot, overlay_polygon

    print(f"Loading samples from {args.sheet} ...")
    samples = load_samples(args.sheet, args.data, subsample=args.subsample)
    print(f"  loaded {len(samples)} samples: {', '.join(s.sample_id for s in samples)}")

    populations = seed_populations(samples)

    # Track per-sample counts through the hierarchy (start at total events).
    counts = {p["sample_id"]: {"n_total": len(p["events"])} for p in populations}
    identity = {p["sample_id"]: {"filename": p["filename"], **p["conditions"]} for p in populations}

    for stage in stages:
        gpath = gate_file(stage)
        if gpath.exists() and not args.redraw:
            gate = load_gate(gpath)
            vertices, x, y = gate["vertices"], gate["x_channel"], gate["y_channel"]
            print(f"[{stage.name}] loaded gate from {gpath} ({len(vertices)} vertices).")
        else:
            x, y = stage.x, stage.y
            print(f"[{stage.name}] pooling {args.per_sample} events/file for {y} vs {x} ...")
            px, py = pool_events(populations, x, y, per_sample=args.per_sample)
            print(f"[{stage.name}] draw the gate (Enter when done) ...")
            vertices = draw_polygon_gate(px, py, x_channel=x, y_channel=y)
            save_gate(vertices, gpath, x, y, name=stage.name)
            print(f"[{stage.name}] saved gate ({len(vertices)} vertices) to {gpath}.")

        # Portable FlowKit gate (for GatingML / FlowJo export).
        build_flowkit_polygon_gate(stage.name, x, y, vertices)

        # Density overlay PNG for this stage (parent population + gate).
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        from cytoflow_vis.style import rc

        ppx, ppy = pool_events(populations, x, y, per_sample=args.per_sample)
        with mpl.rc_context(rc()):
            fig, ax = plt.subplots(figsize=(9, 9))
            density_plot(ppx, ppy, ax=ax, xlabel=x, ylabel=y, title=f"{stage.name} gate",
                         cmap="viridis", colorbar=True)
            overlay_polygon(ax, vertices)
            fig.savefig(out_dir / f"gate_overlay_{stage.name}.png")
            plt.close(fig)

        # Apply: this stage's output becomes the next stage's parent.
        populations = apply_gate(populations, vertices, x, y)
        for p in populations:
            counts[p["sample_id"]][f"n_{stage.name}"] = len(p["events"])

    # Summary table: counts and % of parent through the hierarchy.
    parent_chain = ["total"] + [s.name for s in stages]
    rows = []
    for sid in counts:
        row = {"sample_id": sid, **identity[sid], "n_total": counts[sid]["n_total"]}
        prev = counts[sid]["n_total"]
        for stage in stages:
            n = counts[sid][f"n_{stage.name}"]
            row[f"n_{stage.name}"] = n
            row[f"pct_{stage.name}_of_parent"] = round(100.0 * n / prev, 2) if prev else 0.0
            prev = n
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary_path = out_dir / "gating_summary.csv"
    summary.to_csv(summary_path, index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nWrote summary to {summary_path}.")
    print(f"Wrote {len(stages)} gate overlays to {out_dir}/gate_overlay_*.png.")

    # Final population (singlets) per sample.
    if not args.no_save_events:
        final_name = stages[-1].name
        events_dir = out_dir / f"gated_events_{final_name}"
        events_dir.mkdir(exist_ok=True)
        for p in populations:
            p["events"].to_csv(events_dir / f"{p['sample_id']}.csv", index=False)
        print(f"Wrote {len(populations)} {final_name} event files to {events_dir}/.")


if __name__ == "__main__":
    main()
