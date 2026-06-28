"""`compute-spillover`: build a compensation matrix from single-stain controls.

For legacy FCS files with no embedded $SPILLOVER. Acquire single-stain controls
at the SAME instrument settings as the legacy data, then:

    uv run compute-spillover --controls controls.csv --data controls/ \
        --gates results/ --out spillover.csv

The controls CSV maps each control file to its primary detector:

    filename,channel
    annexin_fitc.fcs,BL1-A
    7aad.fcs,RL1-A

Negative baseline (auto by default):

  * in-tube  -- stain a mix of live + positive cells per tube; the tube's own
    live (negative) cells are the matched baseline. Gold standard; used when no
    unstained control is listed.
  * universal -- add a row with channel ``unstained`` for a separate
    universal-negative control.

Feed the resulting matrix to your analysis config:  compensation = "spillover.csv"
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compute a spillover matrix from single-stain controls.")
    p.add_argument("--controls", required=True, help="Controls CSV (columns: filename, channel).")
    p.add_argument("--data", required=True, help="Directory containing the control FCS files.")
    p.add_argument("--out", default="spillover.csv", help="Output matrix CSV (default: spillover.csv).")
    p.add_argument("--gates", default=None,
                   help="Gate dir; if given, apply cells+singlets gates to the controls (use for stained cells).")
    p.add_argument("--gate-names", default="cells,singlets", help="Ordered gate stage names to replay.")
    p.add_argument("--negative", choices=["auto", "in-tube", "universal"], default="auto",
                   help="Negative baseline: in-tube per-tube split, universal unstained, or auto-detect.")
    p.add_argument("--unstained-label", default="unstained", help="channel value marking the unstained control.")
    p.add_argument("--positive-percentile", type=float, default=99.9,
                   help="Universal-negative percentile defining the positive threshold per detector.")
    p.add_argument("--subsample", type=int, default=50000, help="Events FlowKit stores per control file.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    from flowsmith import fluorescence as fl
    from flowsmith import spillover as sp
    from flowsmith.gating import apply_saved_gates, seed_populations
    from flowsmith.io import load_samples

    print(f"Loading controls from {args.controls} ...")
    samples = load_samples(args.controls, args.data, subsample=args.subsample)
    populations = seed_populations(samples)

    if "channel" not in populations[0]["conditions"]:
        raise SystemExit("controls CSV must have a 'channel' column (primary detector or 'unstained').")

    if args.gates:
        gates_dir = Path(args.gates)
        gate_paths = [gates_dir / f"{n.strip()}_gate.json" for n in args.gate_names.split(",")]
        missing = [str(g) for g in gate_paths if not g.exists()]
        if missing:
            raise SystemExit("Missing gate file(s): " + ", ".join(missing))
        print(f"Gating controls with: {', '.join(g.name for g in gate_paths)} ...")
        populations = apply_saved_gates(populations, gate_paths)

    unstained, single_stains, channels = None, {}, []
    for p in populations:
        ch = p["conditions"]["channel"]
        if ch == args.unstained_label:
            unstained = p["events"]
        else:
            single_stains[ch] = p["events"]
            channels.append(ch)
    if not channels:
        raise SystemExit("no single-stain controls found.")
    available = list(single_stains[channels[0]].columns)
    for ch in channels:
        if ch not in available:
            raise SystemExit(f"detector {ch!r} not found in the FCS channels: {available}")

    # Resolve the negative-baseline method.
    method = args.negative
    if method == "auto":
        method = "universal" if unstained is not None else "in-tube"
    if method == "universal" and unstained is None:
        raise SystemExit(f"--negative universal requires an unstained control (channel == {args.unstained_label!r}).")
    use_unstained = unstained if method == "universal" else None

    extra = f"unstained: {len(unstained)} events" if use_unstained is not None else "in-tube negatives"
    print(f"Detectors: {channels} | negative: {method} ({extra})")
    xform = fl.make_logicle()
    matrix, diag = sp.compute_spillover_matrix(
        single_stains, channels, xform, unstained=use_unstained,
        positive_percentile=args.positive_percentile,
    )

    print("\nSpillover matrix (rows = fluorochrome, cols = detector):")
    print(matrix.round(4).to_string())
    print("\nPer-control diagnostics:")
    for det, d in diag.items():
        neg = f"  negatives={d['n_negative']:>7d}" if d["n_negative"] is not None else ""
        print(f"  {det:8s} positives={d['n_positive']:>7d}{neg}  primary signal={d['primary_separation']}")

    sp.write_matrix_csv(matrix, args.out)
    print(f"\nWrote matrix to {args.out}.  Use it with:  compensation = \"{args.out}\"")


if __name__ == "__main__":
    main()
