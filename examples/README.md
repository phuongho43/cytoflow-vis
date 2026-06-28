# Example configs

Ready-made, annotated `analyze` configs for common flow-cytometry experiments.
Each subfolder is a self-contained template: a `*.toml` config plus an example
`samples.csv` showing the condition columns it expects. Copy a folder, drop your
FCS files into its `data/` directory, adjust the channels/conditions, and run.

Every template assumes the standard two-step flowsmith workflow — draw the
shared FSC/SSC gates once, then run the config-driven analysis:

```bash
cd examples/lentiviral_titration
uv run gate-cells --sheet samples.csv --data data --out results   # once
uv run analyze titration.toml                                     # the analyses
```

Channels in the templates use Attune-style names (`BL1-A` = GFP). Rename them to
your instrument's channels and update the `[channel_labels]` table to match.

## Available templates

| Folder | Experiment | Key readout |
|--------|------------|-------------|
| [`lentiviral_titration/`](lentiviral_titration/) | Serial-dilution titration of a fluorescent-reporter lentivirus | functional titer (TU/mL) + % reporter+ vs virus volume |

More templates (CAR / activation-marker panels, dose-response titrations, …) to
follow.
