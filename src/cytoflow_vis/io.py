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


@dataclass
class LoadedSample:
    """One FCS file plus its experimental conditions."""

    sample_id: str
    filename: str
    conditions: dict
    fcs: fk.Sample
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    def events(self, channels: list[str] | None = None, subsample: bool = False) -> pd.DataFrame:
        """Return events as a DataFrame with flat PnN channel columns.

        ``subsample=False`` returns every event (use this for applying a gate);
        ``subsample=True`` returns FlowKit's stored random subsample (use this
        for fast plotting of a single file).
        """
        key = ("sub" if subsample else "all")
        if key not in self._cache:
            df = self.fcs.as_dataframe(
                source="raw", subsample=subsample, col_multi_index=False
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


def load_samples(
    sheet_path: str | Path,
    data_dir: str | Path,
    subsample: int = 10000,
) -> list[LoadedSample]:
    """Load every FCS file referenced by the sample sheet.

    ``subsample`` sets how many events FlowKit stores as the per-file random
    subsample (used for drawing/plotting); full event data is always available
    via ``LoadedSample.events(subsample=False)``.
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
        conditions = {c: row[c] for c in cond_cols}
        loaded.append(
            LoadedSample(
                sample_id=sample_id,
                filename=str(row[REQUIRED_COLUMN]),
                conditions=conditions,
                fcs=fcs,
            )
        )
    return loaded
