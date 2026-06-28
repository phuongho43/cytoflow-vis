"""Loading FCS files and the experiment sample sheet.

The sample sheet (CSV) is the conditions layer that FlowKit itself does not
impose: it maps each FCS file to its experimental metadata. One column named
``filename`` is required; an optional ``sample`` column gives a short label,
and every other column is treated as an experimental condition.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import flowkit as fk
import pandas as pd

REQUIRED_COLUMN = "filename"
LABEL_COLUMN = "sample"

# Config/argument values that mean "use each FCS file's embedded $SPILLOVER".
EMBEDDED_COMPENSATION = {"acquisition", "embedded", "self", "spillover"}
_SPILL_KEYS = ("spillover", "spill", "comp")


@dataclass
class LoadedSample:
    """One FCS file plus its experimental conditions."""

    sample_id: str
    filename: str
    conditions: dict
    fcs: fk.Sample
    compensated: bool = False
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    def events(self, channels: list[str] | None = None, subsample: bool = False) -> pd.DataFrame:
        """Return events as a DataFrame with flat PnN channel columns.

        Returns compensated values when spillover compensation was applied at
        load time (scatter channels pass through unchanged either way), so
        gating is unaffected while fluorescence uses the corrected values.

        ``subsample=False`` returns every event (use this for applying a gate);
        ``subsample=True`` returns FlowKit's stored random subsample (use this
        for fast plotting of a single file).
        """
        source = "comp" if self.compensated else "raw"
        key = (source, "sub" if subsample else "all")
        if key not in self._cache:
            df = self.fcs.as_dataframe(
                source=source, subsample=subsample, col_multi_index=False
            )
            # Flatten any leftover MultiIndex columns to their PnN label.
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = df.columns.get_level_values(0)
            self._cache[key] = df
        df = self._cache[key]
        if channels is not None:
            return df[list(channels)]
        return df

    @property
    def channels(self) -> list[str]:
        return list(self.events(subsample=True).columns)


def load_sample_sheet(sheet_path: str | Path) -> pd.DataFrame:
    """Read and validate the sample sheet CSV."""
    sheet = pd.read_csv(sheet_path)
    if REQUIRED_COLUMN not in sheet.columns:
        raise ValueError(
            f"Sample sheet {sheet_path!s} must have a '{REQUIRED_COLUMN}' column; "
            f"found columns: {list(sheet.columns)}"
        )
    return sheet


def condition_columns(sheet: pd.DataFrame) -> list[str]:
    """The metadata columns (everything but filename/sample)."""
    return [c for c in sheet.columns if c not in (REQUIRED_COLUMN, LABEL_COLUMN)]


def _resolve_spillover(compensation, fcs: fk.Sample, fcs_path: Path):
    """Return the spillover matrix to apply to ``fcs``.

    ``compensation`` is either one of EMBEDDED_COMPENSATION (use the file's own
    $SPILLOVER keyword) or a path/string to an external matrix that FlowKit's
    ``apply_compensation`` understands (CSV file or FCS spill string).
    """
    if isinstance(compensation, str) and compensation.lower() in EMBEDDED_COMPENSATION:
        for key in _SPILL_KEYS:
            spill = fcs.metadata.get(key)
            if spill:
                return spill
        raise ValueError(
            f"{fcs_path.name}: no embedded spillover matrix "
            f"($SPILLOVER/$SPILL) found in the FCS metadata."
        )
    return compensation  # external CSV path / spill string / Matrix / array


def load_samples(
    sheet_path: str | Path,
    data_dir: str | Path,
    subsample: int = 10000,
    compensation=None,
) -> list[LoadedSample]:
    """Load every FCS file referenced by the sample sheet.

    ``subsample`` sets how many events FlowKit stores as the per-file random
    subsample (used for drawing/plotting); full event data is always available
    via ``LoadedSample.events(subsample=False)``.

    ``compensation`` applies spillover compensation to the fluorescence channels
    at load time. It may be one of ``EMBEDDED_COMPENSATION`` (e.g. "acquisition"
    to use each file's own $SPILLOVER matrix) or an external matrix that
    FlowKit's ``apply_compensation`` accepts (CSV file path or spill string).
    """
    sheet = load_sample_sheet(sheet_path)
    data_dir = Path(data_dir)
    cond_cols = condition_columns(sheet)

    loaded: list[LoadedSample] = []
    for _, row in sheet.iterrows():
        fcs_path = data_dir / row[REQUIRED_COLUMN]
        if not fcs_path.exists():
            raise FileNotFoundError(f"FCS file not found: {fcs_path}")
        sample_id = str(row[LABEL_COLUMN]) if LABEL_COLUMN in sheet.columns else fcs_path.stem
        fcs = fk.Sample(str(fcs_path), sample_id=sample_id, subsample=subsample)
        compensated = False
        if compensation is not None:
            fcs.apply_compensation(_resolve_spillover(compensation, fcs, fcs_path))
            compensated = True
        conditions = {c: row[c] for c in cond_cols}
        loaded.append(
            LoadedSample(
                sample_id=sample_id,
                filename=str(row[REQUIRED_COLUMN]),
                conditions=conditions,
                fcs=fcs,
                compensated=compensated,
            )
        )
    return loaded
