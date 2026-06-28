# Lentiviral titration — protocol & analysis guide

A step-by-step guide to measuring the **functional titer** (transducing units per
mL, TU/mL) of a lentiviral prep carrying a fluorescent reporter (GFP here), and
to running the flowsmith pipeline on the resulting FCS files. Aimed at a
bench biologist; no prior flowsmith experience assumed.

The idea: transduce a **known number of cells** with a **serial dilution** of the
virus, measure the **% reporter-positive** cells by flow a few days later, and
back-calculate the titer from the wells that are dilute enough that most positive
cells carry a single integration.

---

## Part A — Wet lab

> **Biosafety:** lentivirus is handled at **BSL-2**. Follow your institution's
> biosafety approval, use a Class II cabinet, and decontaminate all
> virus-contacted plasticware/media.

**You will need:** your target cell line, the lentiviral prep, a 24-well plate,
polybrene (optional), and a flow cytometer with a laser/filter for your reporter
(GFP → 488 nm excitation, ~510–530 nm emission; "BL1"/"FITC" channel).

1. **Seed a fixed, counted number of cells per well (day 0).** Use the cell type
   you actually care about (titer is **cell-type dependent**). Count carefully —
   e.g. **1 × 10⁵ cells/well** in a 24-well plate. **Write this number down: it is
   `cells_seeded` in the analysis and the titer scales linearly with it.**

2. **Transduce with a serial dilution (day 1).** When cells are ~50–70%
   confluent, add a dilution series of the virus, e.g. **0, 1, 2.5, 5, 10, 25,
   50, 100 µL** per well (top up to a constant total volume with medium). The
   **0 µL well is the uninfected/mock control** — it sets the GFP+ threshold and
   is essential. Run the series in **duplicate or triplicate** (biological
   replicates). Optionally add **polybrene (4–8 µg/mL)** to aid uptake; keep it
   identical across wells.

3. **Incubate 48–72 h** to let the reporter express. Aim to keep the highest
   useful wells **below ~30% positive** if you can — the titer comes from the
   dilute, low-MOI wells, so a few clearly sub-saturating dilutions matter more
   than the saturated ones. Avoid letting cells overgrow.

4. **Harvest and read by flow.** Detach (trypsin/Accutase), wash, resuspend in
   FACS buffer (optionally fix in 1–2% PFA). Acquire **all wells with identical
   cytometer settings** (same PMT voltages/gains). Record ≥10,000 events/well.

5. **Export one FCS file per well.**

---

## Part B — Analysis with flowsmith

This template folder (`examples/lentiviral_titration/`) already contains the two
files you edit: **`samples.csv`** (the well → conditions map) and
**`titration.toml`** (the analysis config). Run the commands **from inside this
folder**.

### B0. One-time setup

```bash
uv sync          # install flowsmith and its dependencies (run once)
```

### B1. Drop in your data and describe the wells

- Put your exported FCS files in a **`data/`** subfolder here.
- Edit **`samples.csv`** so every row points to a real file. Columns:

  | column | meaning |
  |--------|---------|
  | `filename` | the FCS file name in `data/` (required) |
  | `sample` | a short label for the well |
  | `virus_uL` | virus volume added to that well (the dilution series) |
  | `replicate` | biological replicate number (1, 2, …) |

  The provided `samples.csv` already matches the 0–100 µL × 2-replicate series
  above — rename the `filename` entries to your files, or replace the rows.

### B2. Point the config at your reporter and cell count

Open **`titration.toml`** and set two things to match your experiment:

- **`channels = ["BL1-A"]`** → your reporter channel name (e.g. `"FITC-A"`,
  `"FL1-A"`). Update the `[channel_labels]` entry to match so plots read "GFP".
- **`cells_seeded = 100000`** → the number you seeded per well in step A1.

(Leave `control`, `dose`, the linear-range bounds, and `poisson` at their
defaults unless you have a reason to change them — see *Interpreting* below.)

### B3. Draw the gates (interactive, once)

```bash
uv run gate-cells --sheet samples.csv --data data --out results
```

A window opens for each of two gates in turn — first **cells** (SSC-A vs FSC-A,
exclude debris), then **singlets** (FSC-H vs FSC-A, keep the diagonal). For each:
**click to place polygon vertices** around the target population, drag to adjust,
`Esc` to restart, **`Enter`** (or close the window) when happy. The gates are
saved to `results/` as JSON and **reused automatically** on later runs (no
re-clicking; pass `--redraw` to redo them).

### B4. Compute the titer

```bash
uv run analyze titration.toml
```

This replays the gates, applies the logicle transform, and runs the analyses.
The titer line in the console output is your answer, e.g.:

```
titer -> titer.csv, titer.png (mean 5.30e+06 TU/mL over 8 wells)
```

### Outputs (in `results/`)

| File | What it is |
|------|-----------|
| `titer.csv` | Per-well: % positive, the computed `titer_TU_per_mL`, and an `in_range` flag for wells used in the average |
| `titer.png` | Diagnostic plot — per-well titer vs virus volume |
| `fluor_stats.csv` | Per-well MFI + % positive (the underlying numbers) |
| `hist_*.png` | Reporter histogram per virus dose (sanity-check the +/- split) |
| `dose_response.png` | % positive (and MFI) vs virus volume |

---

## Interpreting the result

**How the number is computed.** For each well, the reporter+ fraction `f` is
converted to mean integrations per cell, `MOI = −ln(1 − f)`, and

```
TU/mL = cells_seeded × MOI / (virus volume in mL)
```

The reported titer is the **mean over the in-range wells** (default **5–60%
positive**): below 5% the signal is mostly the control's background, above ~60%
the math gets saturation-sensitive. (Full rationale in
[ASSUMPTIONS.md](../../ASSUMPTIONS.md#viral-titer).)

**Read `titer.png` as a sanity check.** The in-range per-well titers should sit
**flat** across virus volume around the mean line. If they instead **trend with
volume**, the assumptions are strained — see below.

**Troubleshooting:**

| Symptom | Likely cause / fix |
|---------|--------------------|
| Console says *"no wells in 5–60% range"* | All wells too dilute or all saturated. Add lower **or** higher dilutions so some land at 5–60%. |
| Per-well titers **trend down** as volume rises | High wells are saturating; collect more low-dose wells, or lower `linear_max`. |
| Mock (0 µL) well already shows many positives | Background/autofluorescence too high; the threshold is set off the control, so a noisy control inflates everything. Check `hist_*.png`. |
| Titer looks 2–10× off | Re-check `cells_seeded` (titer scales linearly with it) and that `virus_uL`/`volume_unit` are correct. |

**Report it honestly.** State the **cell type**, **`cells_seeded`**, the **time
post-transduction**, that the titer is the mean of the linear-range (single-
integration) wells, and whether polybrene was used — titer is meaningful only
relative to that context.
