"""Compute a spillover (compensation) matrix from single-stain controls.

Spillover of fluorochrome *i* into detector *j* is

    S[i, j] = (pos_median_i[j] - neg_median_i[j])
              / (pos_median_i[i] - neg_median_i[i])

so the diagonal is 1 by construction. Medians are on **raw** (linear) values
because compensation is linear; the logicle transform is only used to separate
positive events from background.

Two ways to obtain the negative baseline:

  * **in-tube** (gold standard): each single-stain tube holds a mix of stained
    (positive) and unstained-equivalent (negative) cells of the same type; the
    positive/negative split is found per tube (Otsu on its primary detector) and
    the tube's own negative is the matched baseline.
  * **universal**: a separate unstained control supplies one negative baseline
    for every detector. Used when an ``unstained`` control is provided.

The result is written as a detector-headed CSV that FlowKit's
``apply_compensation`` reads, so it plugs straight into the ``compensation`` key.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _otsu_threshold(values: np.ndarray, nbins: int = 256) -> float:
    """Otsu's between-class-variance threshold for a bimodal 1D distribution."""
    hist, edges = np.histogram(values, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2.0
    weights = hist.astype(float)
    total = weights.sum()
    if total == 0:
        return float(np.median(values))
    p = weights / total
    omega = np.cumsum(p)
    mu = np.cumsum(p * centers)
    mu_t = mu[-1]
    denom = omega * (1.0 - omega)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b2 = (mu_t * omega - mu) ** 2 / denom
    return float(centers[int(np.nanargmax(sigma_b2))])


def _medians(events: pd.DataFrame, channels: list[str]) -> dict[str, float]:
    return {d: float(np.median(events[d].to_numpy(dtype=float))) for d in channels}


def compute_spillover_matrix(
    single_stains: dict[str, pd.DataFrame],
    channels: list[str],
    xform,
    unstained: pd.DataFrame | None = None,
    positive_percentile: float = 99.9,
    min_events: int = 100,
) -> tuple[pd.DataFrame, dict]:
    """Build an N×N spillover matrix over ``channels``.

    ``single_stains`` maps each primary detector to its (gated) single-stain
    events. If ``unstained`` is given, it is used as the universal negative;
    otherwise an in-tube negative is detected within each single-stain tube.
    Returns the matrix (rows = fluorochromes, columns = detectors, both indexed
    by ``channels``) and a per-control diagnostics dict.
    """
    universal = unstained is not None
    if universal:
        neg_global = _medians(unstained, channels)
        thresholds = {
            d: float(np.percentile(xform.apply(unstained[d].to_numpy(dtype=float)), positive_percentile))
            for d in channels
        }

    rows, diagnostics = {}, {}
    for det in channels:
        if det not in single_stains:
            raise ValueError(f"no single-stain control provided for detector {det!r}")
        events = single_stains[det]
        primary = xform.apply(events[det].to_numpy(dtype=float))

        if universal:
            pos_mask = primary > thresholds[det]
            neg_med = neg_global
            n_neg = None
        else:  # in-tube negative
            thr = _otsu_threshold(primary)
            pos_mask = primary > thr
            neg_events = events[~pos_mask]
            n_neg = int(len(neg_events))
            if n_neg < min_events:
                raise ValueError(
                    f"{det}: only {n_neg} negative events in-tube (< {min_events}); "
                    f"include enough live (unstained) cells in the mix."
                )
            neg_med = _medians(neg_events, channels)

        pos = events[pos_mask]
        n_pos = int(len(pos))
        if n_pos < min_events:
            raise ValueError(
                f"{det}: only {n_pos} positive events (< {min_events}); "
                f"the control is too dim or mislabeled."
            )
        pos_med = _medians(pos, channels)
        denom = pos_med[det] - neg_med[det]
        if denom <= 0:
            raise ValueError(f"{det}: positive median not above negative; cannot compute spillover.")

        rows[det] = {d: (pos_med[d] - neg_med[d]) / denom for d in channels}
        diagnostics[det] = {
            "n_positive": n_pos,
            "n_negative": n_neg,
            "primary_separation": round(denom, 1),
        }

    matrix = pd.DataFrame([rows[d] for d in channels], index=channels, columns=channels)
    return matrix, diagnostics


def write_matrix_csv(matrix: pd.DataFrame, path) -> None:
    """Write the matrix as a detector-headed CSV (FlowKit-readable)."""
    matrix.to_csv(path, index=False)
