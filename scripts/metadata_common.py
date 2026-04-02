from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_NLCD_DIR = REPO_ROOT / "nlcd_data"
DEFAULT_COUNTY_SHP = REPO_ROOT / "us_il_shapefiles" / "cb_2018_us_county_5m.shp"
DEFAULT_CLIMATE_DIR = REPO_ROOT / "climate_data"

DEFAULT_BIOCLIM_FILES = {
    "bio_01": DEFAULT_CLIMATE_DIR / "BIO01_era5_1979-2018_v1.0.nc",
    "bio_04": DEFAULT_CLIMATE_DIR / "BIO04_era5_1979-2018_v1.0.nc",
    "bio_10": DEFAULT_CLIMATE_DIR / "BIO10_era5_1979-2018_v1.0.nc",
    "bio_15": DEFAULT_CLIMATE_DIR / "BIO15_era5_1979-2018_v1.0.nc",
    "bio_17": DEFAULT_CLIMATE_DIR / "BIO17_era5_1979-2018_v1.0.nc",
}


def _pick_column(df: pd.DataFrame, preferred: str | None, candidates: list[str], label: str) -> str:
    if preferred:
        if preferred not in df.columns:
            raise ValueError(f"Requested {label} column '{preferred}' was not found.")
        return preferred

    for column in candidates:
        if column in df.columns:
            return column

    raise ValueError(
        f"Could not find a {label} column. Tried: {', '.join(candidates)}."
    )


def standardize_sample_table(
    df: pd.DataFrame,
    lon_col: str | None = None,
    lat_col: str | None = None,
    year_col: str = "year",
) -> pd.DataFrame:
    normalized = df.copy()

    lon_name = _pick_column(
        normalized,
        lon_col,
        ["decimalLongitude", "long", "longitude", "lon"],
        "longitude",
    )
    lat_name = _pick_column(
        normalized,
        lat_col,
        ["decimalLatitude", "lat", "latitude"],
        "latitude",
    )
    if year_col not in normalized.columns:
        raise ValueError(f"Required year column '{year_col}' was not found.")

    if lon_name != "decimalLongitude":
        normalized["decimalLongitude"] = normalized[lon_name]
    if lat_name != "decimalLatitude":
        normalized["decimalLatitude"] = normalized[lat_name]

    normalized["decimalLongitude"] = pd.to_numeric(
        normalized["decimalLongitude"], errors="coerce"
    )
    normalized["decimalLatitude"] = pd.to_numeric(
        normalized["decimalLatitude"], errors="coerce"
    )
    normalized[year_col] = pd.to_numeric(normalized[year_col], errors="coerce").astype("Int64")

    return normalized


def parse_year_path_overrides(items: Iterable[str] | None) -> dict[int, Path]:
    overrides: dict[int, Path] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(
                f"Invalid override '{item}'. Expected format YEAR=path/to/file.txt."
            )
        year_str, path_str = item.split("=", 1)
        overrides[int(year_str)] = Path(path_str)
    return overrides
