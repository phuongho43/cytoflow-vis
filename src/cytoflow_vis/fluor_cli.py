"""Standalone entry point: fluorescence analysis of the gated singlet population.

Usage (after running `gate-cells` to create the gates):

    uv run analyze-fluor --sheet samples.csv --data data/ --out results/

It reconstructs the singlet population by replaying the saved gates
(``cells_gate.json`` then ``singlets_gate.json``), applies a logicle transform
to the fluorescence channels, and writes:

  * fluor_stats.csv              -- per-sample MFI and % positive
  * fluor_hist_<channel>.png     -- logicle histograms grouped by condition
  * fluor_2d_density.png         -- faceted BL1 vs RL1 density per sample
  * fluor_dose_response.png      -- MFI / % positive vs dose, by condition

This is non-interactive (renders straight to PNG).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Fluorescence analysis of gated singlets.")
    p.add_argument("--sheet", required=True, help="Sample sheet CSV (needs a 'filename' column).")
    p.add_argument("--data", required=True, help="Directory containing the FCS files.")
    p.add_argument("--out", default="results", help="Output directory (default: results/).")
    p.add_argument("--gates", default=None, help="Directory holding the gate JSON files (default: --out).")
    p.add_argument("--gate-names", default="cells,singlets", help="Ordered gate stage names to replay.")
    p.add_argument("--channels", default="BL1-A,RL1-A", help="Fluorescence channels (comma-separated).")
    p.add_argument("--group", default=None, help="Condition column to group/colour by (auto if omitted).")
    p.add_argument("--dose", default=None, help="Numeric condition column for dose-response x-axis (auto).")
    p.add_argument("--control", default=None, help="Sample id used to set the + threshold (auto: lowest dose).")
    p.add_argument("--positive-percentile", type=float, default=99.0, help="Control percentile for + threshold.")
    p.add_argument("--per-sample", type=int, default=20000, help="Events per sample used for plotting.")
    p.add_argument("--subsample", type=int, default=20000, help="Events FlowKit stores per file.")
    p.add_argument("--t", type=float, default=262144.0, help="Logicle param_t (top of raw scale).")
    p.add_argument("--w", type=float, default=0.5, help="Logicle param_w.")
    p.add_argument("--m", type=float, default=4.5, help="Logicle param_m.")
    p.add_argument("--a", type=float, default=0.0, help="Logicle param_a.")
    return p.parse_args(argv)


def _is_numeric(values) -> bool:
    try:
        [float(v) for v in values]
        return True
    except (TypeError, ValueError):
        return False


def detect_columns(populations, group_arg, dose_arg):
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


def dose_determines_group(populations, group_col, dose_col) -> bool:
    """True if each dose value maps to a single group value.

    When the grouping column is redundant with dose (e.g. inducer ``none`` only
    ever occurs at concentration 0 and ``dox`` at the nonzero doses), splitting
    the dose-response by it just orphans the control as its own point. In that
    case we draw one connected curve through every dose instead.
    """
    if group_col is None or dose_col is None:
        return False
    mapping = {}
    for p in populations:
        dose, grp = p["conditions"][dose_col], p["conditions"][group_col]
        if dose in mapping and mapping[dose] != grp:
            return False
        mapping[dose] = grp
    return True


def detect_control(populations, control_arg, dose_col):
    if control_arg is not None:
        return control_arg
    if dose_col is not None:
        return min(populations, key=lambda p: float(p["conditions"][dose_col]))["sample_id"]
    return populations[0]["sample_id"]


def main(argv=None):
    args = parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")  # non-interactive; render straight to PNG
    import matplotlib.pyplot as plt
    import pandas as pd

    from cytoflow_vis import fluorescence as fl
    from cytoflow_vis.gating import apply_saved_gates, seed_populations
    from cytoflow_vis.io import load_samples

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    gates_dir = Path(args.gates) if args.gates else out_dir
    gate_paths = [gates_dir / f"{name.strip()}_gate.json" for name in args.gate_names.split(",")]
    missing = [str(g) for g in gate_paths if not g.exists()]
    if missing:
        raise SystemExit(
            "Missing gate file(s): " + ", ".join(missing) + "\nRun `gate-cells` first to draw them."
        )

    channels = [c.strip() for c in args.channels.split(",")]

    print(f"Loading samples from {args.sheet} ...")
    samples = load_samples(args.sheet, args.data, subsample=args.subsample)
    populations = seed_populations(samples)

    print(f"Replaying gates: {', '.join(g.name for g in gate_paths)} ...")
    singlets = apply_saved_gates(populations, gate_paths)
    for ch in channels:
        if ch not in singlets[0]["events"].columns:
            raise SystemExit(f"Channel {ch!r} not found. Available: {list(singlets[0]['events'].columns)}")
    print(f"  {len(singlets)} samples, e.g. {singlets[0]['sample_id']}: {len(singlets[0]['events'])} singlets")

    group_col, dose_col = detect_columns(singlets, args.group, args.dose)
    control_id = detect_control(singlets, args.control, dose_col)
    print(f"Group column: {group_col} | dose column: {dose_col} | control sample: {control_id}")

    xform = fl.make_logicle(param_t=args.t, param_w=args.w, param_m=args.m, param_a=args.a)

    # 1. Stats: MFI + % positive.
    stats, thresholds = fl.compute_stats(
        singlets, channels, xform, control_id=control_id, positive_percentile=args.positive_percentile
    )
    stats_path = out_dir / "fluor_stats.csv"
    stats.to_csv(stats_path, index=False)
    print("\n" + stats.to_string(index=False))
    print(f"\nWrote stats to {stats_path}.")

    # 2. Histograms by condition (one figure per channel).
    for ch in channels:
        fig, ax = plt.subplots(figsize=(8, 6))
        fl.plot_histograms(
            singlets, ch, xform, ax=ax, group_col=group_col,
            per_sample=args.per_sample, threshold=thresholds.get(ch),
        )
        path = out_dir / f"fluor_hist_{ch.replace('/', '_')}.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {path}.")

    # 3. 2D density (first two channels).
    if len(channels) >= 2:
        fig = fl.plot_2d_density(singlets, channels[0], channels[1], xform, per_sample=args.per_sample)
        path = out_dir / "fluor_2d_density.png"
        fig.savefig(path, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {path}.")

    # 4. Dose-response (MFI and % positive vs dose), if a dose column exists.
    if dose_col is not None:
        # If the group column is just a relabelling of dose, draw one curve
        # through every dose (including the control) instead of splitting it off.
        dr_group = None if dose_determines_group(singlets, group_col, dose_col) else group_col
        metrics = [("MFI", "MFI")] + ([("pct_pos", "% positive")] if thresholds else [])
        nrows = len(metrics)
        fig, axes = plt.subplots(nrows, len(channels), figsize=(6 * len(channels), 5 * nrows), squeeze=False)
        for r, (prefix, ylabel) in enumerate(metrics):
            for c, ch in enumerate(channels):
                ax = axes[r][c]
                fl.plot_dose_response(ax, stats, dose_col, f"{prefix}_{ch}", group_col=dr_group)
                ax.set_ylabel(f"{ylabel} {ch}")
        fig.tight_layout()
        path = out_dir / "fluor_dose_response.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {path}.")
    else:
        print("No numeric dose column detected; skipped dose-response plot.")


if __name__ == "__main__":
    main()
