# cytoflow-vis

Interactive gating and visualization for flow cytometry FCS files, built on
[FlowKit](https://github.com/whitews/FlowKit).

Two commands make up the workflow:

- **`gate-cells`** — draw a **sequential gating hierarchy** of polygon gates and
  apply it across every sample:
  1. **cells** — SSC-A vs FSC-A, excluding low-scatter debris.
  2. **singlets** — FSC-H vs FSC-A, keeping the diagonal and excluding
     doublets/aggregates.

  Each gate is drawn once on a pooled subsample and filtered from the previous
  stage's population (so the singlet gate is drawn on the cells only).

- **`analyze`** — **config-driven** fluorescence analysis of the gated singlets.
  Each experiment is described by a TOML file listing the channels and which
  analysis modules to run (logicle histograms, MFI / % positive, 2D density,
  quadrant, dose-response). See [section 4](#4-fluorescence-analysis).

The gating front-end is shared by every experiment; only the `analyze` config
changes from one experiment to the next.

## Install

```bash
uv sync
```

(Installs FlowKit, matplotlib, and the PySide6 Qt backend for the interactive
window. Requires Python ≥ 3.14.)

## 1. Describe your experiment with a sample sheet

`gate-cells` is driven by a CSV that maps each FCS file to its experimental
conditions. A `filename` column is **required**; an optional `sample` column
gives a short label; every other column is treated as a condition and carried
through to the results.

```csv
filename,sample,genotype,inducer,concentration,replicate
A01.fcs,ctrl,WT,none,0,1
A02.fcs,dox_low,WT,dox,10,1
A03.fcs,dox_high,WT,dox,100,1
```

Put the FCS files in a directory (e.g. `data/`).

## 2. Draw the gates and apply them

```bash
uv run gate-cells --sheet samples.csv --data data/ --out results/
```

For each stage (cells, then singlets):

1. A random subsample of the parent population is **pooled across all files**
   into one density plot, so the gate reflects the whole experiment.
2. An interactive window opens — **click to place polygon vertices** around the
   target population. Drag a vertex to adjust, `Esc` to start over, **`Enter`**
   (or close the window) when done.
3. The gate is applied; its output becomes the next stage's input.

### Outputs (in `--out`)

| File | Contents |
|------|----------|
| `cells_gate.json`, `singlets_gate.json` | The polygon vertices + axes per stage (re-apply without re-clicking) |
| `gating_summary.csv` | Per-sample event counts and **% of parent** at each stage, plus all conditions |
| `gated_events_singlets/<sample>.csv` | The final (singlet) events for each sample |
| `gate_overlay_cells.png`, `gate_overlay_singlets.png` | Density plot + gate overlay per stage |

## 3. Re-apply saved gates (no clicking)

Gates are saved per stage as `<stage>_gate.json` in the output directory (or
`--gates DIR`). On the next run, any stage whose gate file already exists is
reused automatically — only missing gates open a window:

```bash
uv run gate-cells --sheet samples.csv --data data/ --out results/
```

Pass `--redraw` to re-draw every gate from scratch.

### Useful options

| Flag | Default | Purpose |
|------|---------|---------|
| `--gates` | `--out` | Directory holding the `<stage>_gate.json` files |
| `--redraw` | off | Re-draw every gate even if its file exists |
| `--per-sample` | `5000` | Events per file pooled into each draw plot |
| `--subsample` | `20000` | Events FlowKit stores per file |
| `--no-save-events` | off | Skip writing the final per-sample CSVs |

## 4. Fluorescence analysis

The fluorescence stage is **config-driven** so it adapts to each experiment.
You write one TOML file describing the channels and which analyses to run, then:

```bash
uv run analyze experiment.toml
```

It replays the saved gates to reconstruct the singlets, applies a **logicle**
transform, and runs each listed analysis. Paths in the config are relative to
the config file's own location.

```toml
sheet = "samples.csv"
data  = "data"
out   = "results"
channels = ["BL1-A", "RL1-A"]   # default: all non-scatter channels
control  = "light_0"             # negative control for + thresholds (auto: lowest dose)
compensation = "acquisition"     # apply spillover compensation (see below); omit for none
# group / dose columns are auto-detected; override with group = "...", dose = "..."

[[analysis]]
kind = "mfi_pct"                 # MFI + % positive table

[[analysis]]
kind = "histograms"             # one logicle histogram per channel

[[analysis]]
kind = "density_2d"
x = "BL1-A"
y = "RL1-A"

[[analysis]]
kind = "quadrant"               # two thresholds -> four populations
x = "BL1-A"
y = "RL1-A"
labels = { dn = "live", x_pos = "early", dp = "late", y_pos = "necrotic" }

[[analysis]]
kind = "dose_response"          # MFI / % positive vs dose
```

### Built-in analysis kinds

| `kind` | Key params | Outputs |
|--------|------------|---------|
| `mfi_pct` | `channels`, `control`, `positive_percentile` | `fluor_stats.csv` (MFI + % positive) |
| `histograms` | `channels`, `group`, `control` | `hist_<channel>.png` |
| `density_2d` | `x`, `y` | `density_2d.png` (faceted per sample) |
| `quadrant` | `x`, `y`, `labels`, `control`, `dose` | `quadrant.csv`, `*_density.png` (crosshairs + % per quadrant), `*_dose_response.png` (stacked) |
| `dose_response` | `channels`, `dose`, `group`, `control` | `dose_response.png` (MFI + % positive vs dose) |

Every parameter is optional and falls back to the top-level config / auto-detection.
**% positive** and quadrant thresholds are set at `positive_percentile` (default
99th) of the **control** sample's transformed values — use an unstained /
untreated control. Top-level keys: `out`, `gates`, `gate_names`, `subsample`,
`per_sample`, `positive_percentile`, `compensation`, and
`logicle = { param_t = ..., ... }`.

### Spillover compensation

Spillover between fluorochromes is corrected at load time, **before** the
logicle transform and all analyses (scatter channels are untouched, so gating is
unaffected). Set the top-level `compensation` key:

- `compensation = "acquisition"` — use each FCS file's own embedded `$SPILLOVER`
  matrix (the acquisition compensation).
- `compensation = "spillover.csv"` — apply an external matrix file (path relative
  to the config) to every sample.
- omitted — no compensation.

Compensation is essential for panels whose dyes bleed into each other (e.g.
Annexin-V-FITC into the 7-AAD channel): without it, single-positive cells are
misclassified as double-positive, distorting quadrant percentages.

### Adding a new analysis

Write a function in `analysis.py` decorated with `@register("my_kind")` that
takes the `AnalysisContext` plus keyword params and writes to `ctx.out_dir`. It
is then available as `kind = "my_kind"` in any experiment's config.

## Library use

The pieces are importable and compose as a population pipeline:

```python
from cytoflow_vis import load_samples, seed_populations, apply_gate, load_gate

samples = load_samples("samples.csv", "data/")
pops = seed_populations(samples)            # all events, ungated

cells_gate = load_gate("results/cells_gate.json")
cells = apply_gate(pops, cells_gate["vertices"],
                   cells_gate["x_channel"], cells_gate["y_channel"])

singlets_gate = load_gate("results/singlets_gate.json")
singlets = apply_gate(cells, singlets_gate["vertices"],
                      singlets_gate["x_channel"], singlets_gate["y_channel"])

# singlets[i]["events"] is the gated DataFrame for sample i.
```

The same `apply_saved_gates(pops, [".../cells_gate.json", ".../singlets_gate.json"])`
helper replays a gate sequence in one call. Fluorescence analysis then runs on
the resulting populations:

```python
from cytoflow_vis import fluorescence as fl

xform = fl.make_logicle()
stats, thresholds = fl.compute_stats(singlets, ["BL1-A", "RL1-A"], xform, control_id="ctrl")
```

Each gate is also available as a FlowKit `PolygonGate`
(`build_flowkit_polygon_gate`) for export to GatingML / FlowJo.
