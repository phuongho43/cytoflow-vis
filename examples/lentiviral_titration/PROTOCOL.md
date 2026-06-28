# Lentiviral titration — protocol & analysis guide

A step-by-step guide to getting a **functional titer** (transducing units per mL,
TU/mL) for a lentiviral prep carrying a fluorescent reporter (GFP here), and to
running the flowsmith pipeline on the resulting FCS files. Aimed at a bench
biologist; no prior flowsmith experience assumed.

**Scope: in-house QC.** This is a quick, convenient titer to sanity-check a virus
production batch — *not* a publication-grade measurement. We trade rigour for
speed: titer on **HEK293T** (an easy, consistent reporter cell line) rather than
the eventual target cells, and **one well per dilution** (no replicates). The
number is a relative QC readout — good for comparing preps and flagging a failed
production — but the true functional titer on your target cells will differ
(titer is **cell-type dependent**).

The idea: transduce a **known number of cells** with a **serial dilution** of the
virus, measure the **% reporter-positive** cells by flow a few days later, and
back-calculate the titer from the wells dilute enough that most positive cells
carry a single integration.

---

## Part A — Wet lab

> **Biosafety:** lentivirus is handled at **BSL-2**. Follow your institution's
> biosafety approval, use a Class II cabinet, and decontaminate all
> virus-contacted plasticware/media.

**You will need:** HEK293T cells, the lentiviral prep (~100 µL), a 12-well plate,
polybrene (optional), and a flow cytometer with a laser/filter for your reporter
(GFP → 488 nm excitation, ~510–530 nm emission; "BL1"/"FITC" channel).

1. **Seed and transduce in one step.** Plate **1 × 10⁵ HEK293T in 0.5 mL DMEM per
   well** of a 12-well plate and **add the virus to each well at the same time**
   (no separate next-day transduction). Because virus goes on at seeding, **100K
   is the cell count at the moment of transduction** — that is exactly the number
   the titer formula needs, so use `cells_seeded = 100000` and do **not** re-count
   at harvest (see the note below). Optionally add **polybrene (4–8 µg/mL)**,
   identical across wells.

2. **Make the dilution series.** From your 100 µL prep, take **~50 µL** and make a
   **~3-fold (third-log) serial dilution**, then add a fixed volume to each well to
   give **original-prep-equivalent volumes** of **0.01, 0.03, 0.1, 0.3, 1, 3, 10,
   30, 100 µL**, plus a **0 µL mock** — **one well per dilution** (10 wells, fits a
   12-well plate). This nine-step span is the recommended default: it guarantees
   **≥2 wells land in the usable 5–60% positive range for any titer from 10⁶ to
   10⁹ TU/mL**, so you don't need to know the titer beforehand (the analysis keeps
   the in-range wells and ignores the rest). Sub-µL equivalents are made by
   diluting the prep — e.g. 0.01 µL-equivalent = 1 µL of a 1:100 dilution; 10 µL of
   a 1:10 dilution → `virus_uL = 1`. Record the prep-equivalent volume as
   `virus_uL` in `samples.csv`, so the reported titer refers to your **undiluted**
   prep. The **0 µL mock** sets the GFP+ threshold and is essential.

3. **Incubate 48–72 h** to let the reporter express. Aim to keep the highest
   useful wells **below ~30% positive** if you can — the titer comes from the
   dilute, low-MOI wells, so a few clearly sub-saturating dilutions matter more
   than the saturated ones. Avoid letting the fast-growing HEK293T overgrow.

4. **Harvest and read by flow.** Detach (trypsin/Accutase), wash, resuspend in
   FACS buffer (optionally fix in 1–2% PFA). Acquire **all wells with identical
   cytometer settings** (same PMT voltages/gains). Record ≥10,000 events/well.

5. **Export one FCS file per well.**

> **Why `cells_seeded` is the transduction-time count, not the harvest count.**
> The titer is `cells_seeded × MOI / volume`, and `MOI = −ln(1 − f)` is defined
> *per cell present when the virus infects*. The measured % positive is a ratio
> that is preserved as cells divide (a transduced cell's progeny stay GFP+), so
> all the absolute-titer information sits in the cell number at transduction.
> Counting at harvest would multiply in the 48–72 h of HEK293T divisions and
> overestimate the titer several-fold. Using 100K is both correct and more
> reproducible than a harvest count.

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

  The provided `samples.csv` already matches the third-log 0.01–100 µL series
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
titer -> titer.csv, titer.png (mean 5.30e+06 TU/mL over 4 wells)
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

**Treat it as a relative QC number.** This is a HEK293T titer from single wells —
use it to compare preps, track production consistency, and catch a failed batch,
not as an absolute titer for your target cells. When you record it in the lab
notebook, note the **cell type (HEK293T)**, **`cells_seeded`**, **time
post-transduction**, and whether polybrene was used, so batches stay comparable.
If you ever need a publication-grade or target-cell titer, switch to the actual
target cells and run replicates.
