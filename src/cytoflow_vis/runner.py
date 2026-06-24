"""Config-driven analysis runner: the `analyze` entry point.

Reads a per-experiment TOML config describing the data, channels, and the list
of analyses to run, then:

  1. loads the samples and replays the saved gates (cells -> singlets),
  2. builds a shared AnalysisContext (logicle transform + auto-detected
     dose/group/control columns, overridable in the config),
  3. dispatches each ``[[analysis]]`` entry to its registered module.

Example config:

    sheet = "samples.csv"
    data  = "data/"
    out   = "results/"
    channels = ["BL1-A", "RL1-A"]
    control  = "light_0"

    [[analysis]]
    kind = "mfi_pct"

    [[analysis]]
    kind = "quadrant"
    x = "BL1-A"
    y = "RL1-A"
    labels = { dn = "live", x_pos = "early", dp = "late", y_pos = "necrotic" }

Run with:  uv run analyze experiment.toml
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path


def load_config(path: str | Path) -> dict:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Run a config-driven fluorescence analysis.")
    p.add_argument("config", help="Path to the experiment TOML config file.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = load_config(args.config)
    cfg_dir = Path(args.config).resolve().parent

    def resolve(p: str) -> Path:
        # Paths in the config are relative to the config file's location.
        p = Path(p)
        return p if p.is_absolute() else cfg_dir / p

    import matplotlib

    matplotlib.use("Agg")  # non-interactive; render straight to PNG

    from cytoflow_vis import analysis as an
    from cytoflow_vis import fluorescence as fl
    from cytoflow_vis.gating import apply_saved_gates, seed_populations
    from cytoflow_vis.io import load_samples

    if "sheet" not in cfg or "data" not in cfg:
        raise SystemExit("config must define 'sheet' and 'data'.")

    out_dir = resolve(cfg.get("out", "results"))
    out_dir.mkdir(parents=True, exist_ok=True)
    gates_dir = resolve(cfg["gates"]) if "gates" in cfg else out_dir
    gate_names = cfg.get("gate_names", ["cells", "singlets"])
    gate_paths = [gates_dir / f"{name}_gate.json" for name in gate_names]
    missing = [str(g) for g in gate_paths if not g.exists()]
    if missing:
        raise SystemExit("Missing gate file(s): " + ", ".join(missing) + "\nRun `gate-cells` first.")

    channels = cfg.get("channels")
    subsample = int(cfg.get("subsample", 20000))

    # Spillover compensation: a keyword (e.g. "acquisition") uses each file's
    # embedded $SPILLOVER; anything else is an external matrix file (resolved
    # relative to the config).
    from cytoflow_vis.io import EMBEDDED_COMPENSATION

    comp = cfg.get("compensation")
    if comp is not None and str(comp).lower() not in EMBEDDED_COMPENSATION:
        comp = str(resolve(comp))

    print(f"Loading samples from {cfg['sheet']} ...")
    if comp is not None:
        print(f"Applying spillover compensation: {comp}")
    samples = load_samples(
        resolve(cfg["sheet"]), resolve(cfg["data"]), subsample=subsample, compensation=comp
    )
    populations = seed_populations(samples)
    print(f"Replaying gates: {', '.join(g.name for g in gate_paths)} ...")
    singlets = apply_saved_gates(populations, gate_paths)

    if channels is None:  # default: all non-scatter channels
        scatter = {"FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W", "Time"}
        channels = [c for c in singlets[0]["events"].columns if c not in scatter]
    for ch in channels:
        if ch not in singlets[0]["events"].columns:
            raise SystemExit(f"Channel {ch!r} not found. Available: {list(singlets[0]['events'].columns)}")

    group_col, dose_col = an.detect_columns(singlets, cfg.get("group"), cfg.get("dose"))
    control_id = an.detect_control(singlets, cfg.get("control"), dose_col)
    print(f"Channels: {channels} | group: {group_col} | dose: {dose_col} | control: {control_id}")

    ctx = an.AnalysisContext(
        populations=singlets,
        out_dir=out_dir,
        xform=fl.make_logicle(**cfg.get("logicle", {})),
        channels=channels,
        control_id=control_id,
        group_col=group_col,
        dose_col=dose_col,
        per_sample=int(cfg.get("per_sample", 20000)),
        positive_percentile=float(cfg.get("positive_percentile", 99.0)),
    )

    analyses = cfg.get("analysis", [])
    if not analyses:
        print("No [[analysis]] entries in config; nothing to do.")
        return
    print(f"Running {len(analyses)} analyses into {out_dir}/ ...")
    for entry in analyses:
        params = dict(entry)
        kind = params.pop("kind", None)
        if kind is None:
            raise SystemExit("each [[analysis]] needs a 'kind'.")
        print("  " + an.run_analysis(kind, ctx, params))


if __name__ == "__main__":
    main()
