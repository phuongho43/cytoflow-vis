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

- **`analyze-fluor`** — fluorescence analysis of the gated singlets (logicle
  transform, histograms by condition, MFI / % positive, 2D density,
  dose-response). See [section 4](#4-fluorescence-analysis).

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

Once the gates exist, analyse the fluorescence of the singlet population:

```bash
uv run analyze-fluor --sheet samples.csv --data data/ --out results/
```

It replays the saved gates to reconstruct the singlets, applies a **logicle**
transform to the fluorescence channels (default `BL1-A,RL1-A`), and writes:

| File | Contents |
|------|----------|
| `fluor_stats.csv` | Per-sample **MFI** (median, raw scale) and **% positive** per channel, with conditions |
| `fluor_hist_<channel>.png` | Logicle histograms, one line per sample, coloured by condition |
| `fluor_2d_density.png` | Faceted BL1-A vs RL1-A logicle density, one panel per sample |
| `fluor_dose_response.png` | MFI and % positive vs dose, grouped by condition |

**% positive** is computed against a threshold set at the
`--positive-percentile` (default 99th) of a **control** sample's transformed
values — so use an unstained / uninduced control. The dose and grouping columns
are auto-detected from the sample sheet; override with `--dose`, `--group`, and
`--control`.

### Useful options

| Flag | Default | Purpose |
|------|---------|---------|
| `--channels` | `BL1-A,RL1-A` | Fluorescence channels to analyse |
| `--control` | lowest dose | Sample used to set the + threshold |
| `--dose` / `--group` | auto | Condition columns for the x-axis / colour |
| `--positive-percentile` | `99.0` | Control percentile defining "positive" |
| `--t` `--w` `--m` `--a` | logicle std | Logicle transform parameters |

> Note: this step assumes any spillover **compensation** was already applied (or
> isn't needed). Compensation from the FCS `$SPILLOVER` matrix is a possible
> future addition.

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
