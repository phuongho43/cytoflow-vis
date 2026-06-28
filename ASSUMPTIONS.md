# Analysis assumptions & caveats

A living record of the assumptions baked into `flowsmith`'s analyses, so that
outputs are interpreted correctly and reported honestly in figure legends /
methods. Add to this file whenever a new analysis introduces an assumption.

## Statistics

### Error bars

- **Default is mean ± SD** (standard deviation of the replicates) — a
  *descriptive* measure of spread that does not shrink with `n` and never
  overstates precision. Configurable per analysis via `errorbar`:
  - `sd` — standard deviation (default).
  - `sem` — standard error of the mean (`SD/√n`); precision of the mean.
  - `ci` — a **proper 95 % confidence interval using the t-distribution**
    (`t(0.975, n−1) × SEM`), **not** the `1.96 × SEM` normal approximation
    (which understates the interval at small `n`: the correct multiplier is
    ≈4.30 at n=3, ≈2.78 at n=5, → 1.96 only as n→∞).
- A confidence interval describes the **precision of the mean**, not the spread
  of the data; its frequentist meaning is the long-run capture rate of the
  procedure, not "95 % probability the mean is in this interval".
- **Individual replicate points are always plotted** alongside the bar. At the
  small `n` typical of these experiments this matters more than the choice of
  error bar, and is expected by Cell/Nature/Science.

### Significance brackets

- Each bracket is an **independent two-sample Welch's t-test**
  (`scipy.stats.ttest_ind(..., equal_var=False)`), drawn as `* ≤ 0.05`,
  `** ≤ 0.01`, `*** ≤ 0.001`, else `ns`.
- **No multiple-comparison correction is applied.** This is fine for a few
  *planned* comparisons; if you draw many brackets, apply a correction
  (Holm/Bonferroni) or use a one-way **ANOVA + post-hoc** test (Tukey/Dunnett)
  instead, and state it. (*Post-hoc* = the follow-up pairwise tests run after an
  ANOVA shows some difference, with the multiple-comparison correction built in.)
- The t-test assumes approximately **normally distributed** values, which cannot
  be verified at small `n` (e.g. n=3) — treat it as an assumption.
- Brackets and bars are independent: the **bar shows spread/precision, the
  bracket shows the test result**. Do not infer significance from error-bar
  overlap (overlapping 95 % CIs can still differ at p<0.05).

### What to put in the figure legend

Error-bar type, **n**, **biological vs technical** replicates, the **test** used
(one/two-tailed, paired/unpaired) and any **multiple-comparison correction**.

### Replicates

- Display plots (histograms, density/quadrant facets) show **one representative
  replicate** per condition — the one whose median signal is closest to the
  group median (`representative_population`); they are illustrative, not the
  statistic. Quantitative plots (dose-response, categorical, quadrant
  dose-response) **aggregate across replicates** (mean ± error over the
  individual points). No pooling of biological replicates for stats
  (no pseudoreplication); the sample-sheet `replicate` column defines them.

## Gating, transforms & thresholds

- **Logicle (biexponential) transform** is applied to fluorescence channels for
  display and thresholding (defaults `param_t=262144, param_w=0.5, param_m=4.5,
  param_a=0`, overridable via the config `logicle` table). Scatter channels are
  left linear.
- **% positive and quadrant thresholds** are set at `positive_percentile`
  (default 99th) of a **control** sample's transformed values. The control
  **must be a true negative** (unstained / untreated / lowest dose) — if it
  already contains positive events, the threshold lands inside a positive
  population and mis-classifies everything downstream.
- **MFI** is the median of the **raw** (untransformed) values — the standard
  reporting convention.
- **Subsampling**: each sample is subsampled to `per_sample` events. This also
  makes the density `events / bin` scale comparable across panels (equal N).

## Density & quadrant facets

- Faceted density is a **raw 2D histogram** (the cytometry-standard
  pseudocolour plot) with a single shared **absolute `events / bin`** colour
  scale across panels; comparability relies on the equal-N subsampling above.
- The **quadrant** gate is **binary** (above/below each threshold). It reports
  *what fraction crossed* a threshold, not *how strongly* — two conditions can
  both read "100 % positive" while differing greatly in MFI. Read the quadrant
  % and the MFI together.

## Dose-response

- The dose axis is **symlog** with the linear-region threshold set to the
  smallest positive dose, so a zero dose has a place and fractional doses spread
  out instead of bunching at 0.
- The connecting line joins the **means** — it is **not a fitted model**. There
  is currently **no curve fitting / EC50** (a 4-parameter-logistic fit is on the
  [roadmap](ROADMAP.md)); do not read an EC50 off the connected points.

## Viral titer

- **Titer is a derived quantity, not a measurement.**
  `TU/mL = cells_seeded × transduced_per_cell / volume_mL`. `cells_seeded` (target
  cells per well *at the time of transduction*) is supplied by the user — flow
  cytometry cannot recover it — so the titer is only as accurate as that number.
- **Poisson single-hit correction (default).** Each well's reporter+ fraction `f`
  is converted to the mean integrations per cell before scaling,
  `MOI = −ln(1 − f)`, undoing the bias from cells taking up several particles. It
  assumes integrations are **independent and Poisson-distributed**, cells are
  **equally infectable**, and one integration suffices to score a cell positive —
  approximations that hold best at low-to-moderate MOI. `poisson = false` uses the
  raw fraction `f`, which is unbiased *only* at low MOI and systematically
  **underestimates** the titer as wells saturate.
- **Only the linear-range wells are averaged** (`linear_min`–`linear_max`,
  default **5–60 %** positive). The **lower** bound drops wells dominated by the
  control's false-positive background: the positive threshold sits at a high
  percentile of the control (≈1 % false positives by construction), and at small
  `f` that background is a large *relative* error that `−ln(1 − f)` then
  amplifies. The **upper** bound stays clear of the ≳90 % saturation regime, where
  `−ln(1 − f)` is steep (`d/df = 1/(1−f)`) and hypersensitive to small
  fraction-measurement errors (and where receptor saturation / toxicity break the
  Poisson assumptions). The reported titer is the **mean** over the in-range wells.
- **Read `titer.png` as a diagnostic.** In a clean titration the per-well titers
  sit **flat** across virus volume; a trend with volume means the window or the
  model assumptions are off (raw mode in particular trends downward as wells
  saturate). Reported value = the dashed mean line.

---

*Add new assumptions here as analyses are added (e.g. proliferation modelling,
curve fitting, clustering). See also [ROADMAP.md](ROADMAP.md).*
