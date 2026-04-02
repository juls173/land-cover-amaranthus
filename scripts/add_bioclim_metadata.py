from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import xarray as xr

from metadata_common import DEFAULT_BIOCLIM_FILES, standardize_sample_table


def infer_variable_name(ds: xr.Dataset, output_column: str) -> str:
    stem_token = output_column.replace("_", "").upper()
    candidates = [
        stem_token,
        output_column.upper(),
        output_column.lower(),
        output_column,
    ]
    for candidate in candidates:
        if candidate in ds.data_vars:
            return candidate

    data_vars = list(ds.data_vars)
    if len(data_vars) != 1:
        raise KeyError(
            f"Could not infer climate variable for '{output_column}'. Available variables: {data_vars}"
        )
    return data_vars[0]


def infer_lat_lon_names(ds: xr.Dataset) -> tuple[str, str]:
    lat_name = "lat" if "lat" in ds.coords else "latitude" if "latitude" in ds.coords else None
    lon_name = "lon" if "lon" in ds.coords else "longitude" if "longitude" in ds.coords else None
    if lat_name is None or lon_name is None:
        raise KeyError(f"Could not find lat/lon coordinates in dataset. Found: {list(ds.coords)}")
    return lat_name, lon_name


def select_year_slice(
    data_array: xr.DataArray,
    years: pd.Series,
    year_selection: str = "nearest",
) -> xr.DataArray:
    years_da = xr.DataArray(years.astype(int).values, dims="rows")
    if "year" in data_array.coords:
        if year_selection == "nearest":
            return data_array.sel(year=years_da, method="nearest")
        return data_array.sel(year=years_da)

    if "time" in data_array.coords or "time" in data_array.dims:
        annual = data_array.groupby("time.year").mean("time")
        if year_selection == "nearest":
            return annual.sel(year=years_da, method="nearest")
        return annual.sel(year=years_da)

    raise KeyError("No 'year' or 'time' coordinate found in climate dataset.")


def sample_bioclim_from_dataset(
    df: pd.DataFrame,
    output_column: str,
    nc_path: str | Path,
    year_col: str = "year",
    interp_method: str = "nearest",
    year_selection: str = "nearest",
    engine: str | None = None,
) -> pd.Series:
    with xr.open_dataset(nc_path, engine=engine) as ds:
        variable_name = infer_variable_name(ds, output_column)
        data_array = select_year_slice(ds[variable_name], df[year_col], year_selection=year_selection)
        lat_name, lon_name = infer_lat_lon_names(ds)

        longitudes = df["decimalLongitude"].copy()
        if float(ds[lon_name].max()) > 180 and float(longitudes.min(skipna=True)) < 0:
            longitudes = (longitudes + 360) % 360

        lats = xr.DataArray(df["decimalLatitude"].values, dims="rows")
        lons = xr.DataArray(longitudes.values, dims="rows")

        if interp_method == "nearest":
            values = data_array.sel({lat_name: lats, lon_name: lons}, method="nearest")
        else:
            values = data_array.interp({lat_name: lats, lon_name: lons}, method=interp_method)

        return pd.Series(values.values, index=df.index, name=output_column)


def add_bioclim_metadata(
    df: pd.DataFrame,
    variable_to_path: dict[str, str | Path] | None = None,
    year_col: str = "year",
    lon_col: str | None = None,
    lat_col: str | None = None,
    interp_method: str = "nearest",
    year_selection: str = "nearest",
    engine: str | None = None,
) -> pd.DataFrame:
    result = standardize_sample_table(df, lon_col=lon_col, lat_col=lat_col, year_col=year_col)
    variable_to_path = variable_to_path or DEFAULT_BIOCLIM_FILES

    valid_mask = (
        result["decimalLongitude"].notna()
        & result["decimalLatitude"].notna()
        & result[year_col].notna()
    )
    valid = result.loc[valid_mask].copy()

    for output_column, nc_path in variable_to_path.items():
        result[output_column] = pd.NA
        sampled = sample_bioclim_from_dataset(
            valid,
            output_column=output_column,
            nc_path=nc_path,
            year_col=year_col,
            interp_method=interp_method,
            year_selection=year_selection,
            engine=engine,
        )
        result.loc[valid.index, output_column] = sampled

    return result


def parse_variable_path_overrides(items: list[str] | None) -> dict[str, Path]:
    overrides: dict[str, Path] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(
                f"Invalid climate override '{item}'. Expected format bio_15=path/to/file.nc."
            )
        variable, path = item.split("=", 1)
        overrides[variable.strip()] = Path(path)
    return overrides


def main() -> None:
    parser = argparse.ArgumentParser(description="Add bioclim variables from NetCDF rasters to a sample table.")
    parser.add_argument("--input", "-i", required=True, help="Input CSV with coordinates and year.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path.")
    parser.add_argument(
        "--variable",
        action="append",
        default=[],
        metavar="BIO_COL=PATH",
        help="Override or add a climate source, e.g. bio_15=climate_data/BIO15_...nc",
    )
    parser.add_argument("--year-col", default="year", help="Name of the year column.")
    parser.add_argument("--lon-col", default=None, help="Override longitude column name.")
    parser.add_argument("--lat-col", default=None, help="Override latitude column name.")
    parser.add_argument(
        "--interp-method",
        default="nearest",
        help="Spatial interpolation method passed to xarray (default: nearest).",
    )
    parser.add_argument(
        "--year-selection",
        choices=["nearest", "exact"],
        default="nearest",
        help="How to select years from the climate dataset (default: nearest).",
    )
    parser.add_argument("--engine", default=None, help="Optional xarray backend engine.")
    args = parser.parse_args()

    variable_to_path = dict(DEFAULT_BIOCLIM_FILES)
    variable_to_path.update(parse_variable_path_overrides(args.variable))

    df = pd.read_csv(args.input)
    enriched = add_bioclim_metadata(
        df,
        variable_to_path=variable_to_path,
        year_col=args.year_col,
        lon_col=args.lon_col,
        lat_col=args.lat_col,
        interp_method=args.interp_method,
        year_selection=args.year_selection,
        engine=args.engine,
    )
    enriched.to_csv(args.output, index=False)
    print(f"Wrote bioclim metadata to: {args.output}")


if __name__ == "__main__":
    main()
