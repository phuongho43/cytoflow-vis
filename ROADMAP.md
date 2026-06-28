# Roadmap (tentative)

A loose, prioritized list of analyses and visualizations that are common in
flow-cytometry figures (incl. Cell/Nature/Science) and that `cytoflow-vis` does
not yet cover. Priorities are weighted toward this project's scope — conventional
/ few-colour immunoengineering flow (CAR-T, reporters, titrations) — rather than
high-parameter spectral or CyTOF panels. Nothing here is committed work; it's a
direction sketch and ordering is open to change.

## Current capabilities (for context)

- Sequential gating (cells → singlets), polygon gates saved/replayed, gate-overlay plots
- Spillover compensation from single-stain controls
- Logicle (biexponential) display transform
- 1D distributions as **ridgelines** (per-condition, threshold + % positive)
- 2D density facets (raw hist2d, shared `events/bin` scale, **row/col grid layout**)
- Quadrant analysis (two thresholds → 4 populations) + faceted density + quadrant dose-response
- **Dose-response** (MFI / % positive vs dose, multi-group lines)
- **Categorical comparison** (unordered groups, subgroup dodge, mean ± 95% CI, t-test significance brackets)
- Per-sample MFI / % positive tables

## Near-term, high value (fit current scope)

- [ ] **Dose-response curve fitting (4-parameter logistic / Hill) + EC50/IC50.**
      We already plot and connect dose-response points; top papers *fit* the
      sigmoid and report EC50 ± CI. Natural extension of `dose_response`.
- [ ] **Proliferation / dye-dilution (CFSE, CTV).** Model division peaks →
      % divided, division index, proliferation index. Core CAR-T / T-cell
      functional readout.
- [ ] **Population-frequency composition bars.** "% of each population across
      conditions." Two forms (see notes): a *stacked* composition view and a
      *grouped* per-population view. The grouped view is ~free (it is
      `plot_categorical` fed frequencies); the stacked bar is the only new plot.
      Needs an upstream step to tabulate per-population frequencies (we have
      `quadrant_stats` for the 4-quadrant case; generalize to arbitrary gates).
- [ ] **Cytotoxicity / functional dose-response.** % specific lysis vs E:T ratio
      — a killing-assay variant of `dose_response`. Common in CAR-T papers.
- [ ] **Gating hierarchy / N-level sequential gating + frequency-of-parent.**
      Arbitrary-depth gating (live → lineage → subset) with a population tree and
      % at each node. We currently stop at cells → singlets.
- [ ] **Violin / box option for the categorical comparison.** Alternative to the
      strip + CI when a fuller distribution view is wanted.
- [ ] **Significance testing: ANOVA + post-hoc, with auto-selected tests.**
      Current significance brackets are *uncorrected pairwise Welch t-tests* —
      fine for two groups or a few planned contrasts, but the wrong tool once
      there are many comparisons or a two-factor design. Add:
      - one-way **ANOVA** (≥3 groups) → post-hoc (**Tukey** all-pairwise, or
        **Dunnett** vs a named control);
      - two-way **ANOVA** for `group × subgroup` designs, reporting the main
        effects *and the interaction* (the interaction is often the real
        question, e.g. "does CAR-A respond to stimulation more than CAR-B?");
      - multiple-comparison **correction** (Holm/Bonferroni) for sets of
        planned pairwise brackets.

      **Auto-selection** is feasible from structure the package already tracks
      (group/subgroup columns, their data types, level counts, whether a control
      is named): 1 factor/2 levels → t-test; 1 factor/≥3 → one-way ANOVA; two
      factors → two-way ANOVA + interaction; explicit `sig` pairs → planned
      comparisons + correction. Follow the existing "auto-detect a sensible
      default, but report it and allow override" pattern. **Do not** auto-switch
      parametric↔non-parametric from underpowered small-n normality tests
      (default to robust Welch variants; let the user opt into Mann–Whitney /
      Kruskal–Wallis), **do not** auto-decide paired vs unpaired (it's metadata
      the package can't infer — default unpaired, make it explicit), and treat a
      numeric factor as categorical-vs-trend a user choice. See
      [ASSUMPTIONS.md](ASSUMPTIONS.md).

## Considered but deprioritized

- **Overlaid histograms** (multiple conditions superimposed on one axis).
  Ridgelines are superior for the many-group series that dominates this package,
  and our ridgelines already carry the threshold line + % positive annotation
  that the 2–3 group overlay is usually drafted for. Real value is limited to the
  narrow few-group *direct-overlap* contrast (e.g. stained vs FMO). Low priority.

## Higher-parameter analyses (larger efforts; mainly for high-dimensional panels)

- [ ] **Dimensionality reduction** — UMAP / t-SNE embeddings coloured by marker
      or cluster.
- [ ] **Unsupervised clustering** — FlowSOM / Leiden / Phenograph, with
      cluster × marker heatmaps and cluster-abundance comparison (diffcyt-style).
- [ ] **Marker expression heatmaps** — z-scored MFI, sample × marker / population.

These are ubiquitous in modern immunophenotyping papers but add substantial
machinery (embedding/clustering deps) and may be beyond a few-colour workflow.

## Specialized / assay-specific

- [ ] **Cell-cycle analysis** (DNA content → G0/G1, S, G2/M fractions)
- [ ] **Calcium flux / kinetics** (time vs ratio; we have a Time channel but no kinetic plot)
- [ ] **Backgating** (show a gated population's location in parent plots)
- [ ] **Volcano plots** for differential abundance across many markers/clusters
- [ ] **MFI fold-change / log2-normalized output** (normalize to a reference condition)

## Developer experience

- [ ] **Example config templates.** A folder of ready-made TOML configs for common
      experiment types (lentiviral titration, CAR/CD69 activation, dox/Tet
      induction, gating-only, …) so users copy a template instead of writing one
      from scratch. Should showcase `channel_labels`, `group`/`subgroup`/`order`,
      `sig` brackets, the `row`/`col` facet grid, `show_n`, `cmap`, etc.
