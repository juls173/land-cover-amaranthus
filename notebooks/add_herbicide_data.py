# pip install pandas geopandas shapely pyarrow

import re
import pathlib
import pandas as pd
import geopandas as gpd
import sys

from shapely.geometry import Point

# ------------------------------------------------------------------------------------
# CONFIG — update these 3 paths to match your machine
# ------------------------------------------------------------------------------------
SITES_CSV   = r"C:\\Users\\julia\\Downloads\\sites_with_pesticides - Sheet1 (1).csv"
COUNTY_SHP  = r"C:\\Users\\julia\\Documents\\GitHub\\genomics_lulc\\us_il_shapefiles\\cb_2018_us_county_5m.shp"   # any TIGER/Line county shapefile works
EPEST_FILE = r"..\EPest_county_estimates_2018.txt"
EPEST_BASE = "https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/PesticideUseEstimates"
EPEST_URL  = EPEST_BASE + "/EPest.county.estimates.{year}.txt"
EPEST_DIR   = r".\epest_cache"

# Output
OUT_PATH    = r".\sites_with_action_sites.csv"

# Metric choice: "high", "low", or "mean" of EPEST_LOW_KG and EPEST_HIGH_KG
OUT_PATH         = r".\sites_with_action_sites.csv"

# Which EPest estimate to use: "high" | "low" | "mean" (avg of low/high)
VALUE_KIND       = "high"

# Years you want to force to local files. Add more entries if needed.
LOCAL_BY_YEAR = {
    2018: r"..\EPest_county_estimates_2018.txt",     # ← you said only 2018 is local; add others here if needed
}

# Toggle simple progress logging
VERBOSE = True

# ---- Paste your full lookup dict here (exactly as you shared earlier) ----
epest_herbicide_moa = {
    # ALS inhibitors (Group 2) - 35 compounds
    'BENSULFURON': 'ALS inhibitor (Group 2)',
    'CHLORIMURON': 'ALS inhibitor (Group 2)',
    'CHLORSULFURON': 'ALS inhibitor (Group 2)',
    'CLORANSULAM-METHYL': 'ALS inhibitor (Group 2)',
    'DICLOSULAM': 'ALS inhibitor (Group 2)',
    'FLAZASULFURON': 'ALS inhibitor (Group 2)',
    'FLORASULAM': 'ALS inhibitor (Group 2)',
    'FLUCARBAZONE': 'ALS inhibitor (Group 2)',
    'FORAMSULFURON': 'ALS inhibitor (Group 2)',
    'HALOSULFURON': 'ALS inhibitor (Group 2)',
    'IMAZAMOX': 'ALS inhibitor (Group 2)',
    'IMAZAPIC': 'ALS inhibitor (Group 2)',
    'IMAZAPYR': 'ALS inhibitor (Group 2)',
    'IMAZAQUIN': 'ALS inhibitor (Group 2)',
    'IMAZETHAPYR': 'ALS inhibitor (Group 2)',
    'IMAZOSULFURON': 'ALS inhibitor (Group 2)',
    'IODOSULFURON': 'ALS inhibitor (Group 2)',
    'MESOSULFURON': 'ALS inhibitor (Group 2)',
    'METSULFURON': 'ALS inhibitor (Group 2)',
    'NICOSULFURON': 'ALS inhibitor (Group 2)',
    'ORTHOSULFAMURON': 'ALS inhibitor (Group 2)',
    'PENOXSULAM': 'ALS inhibitor (Group 2)',
    'PRIMISULFURON': 'ALS inhibitor (Group 2)',
    'PROPOXYCARBAZONE': 'ALS inhibitor (Group 2)',
    'PROSULFURON': 'ALS inhibitor (Group 2)',
    'PYRITHIOBAC-SODIUM': 'ALS inhibitor (Group 2)',
    'PYROXSULAM': 'ALS inhibitor (Group 2)',
    'RIMSULFURON': 'ALS inhibitor (Group 2)',
    'SULFOSULFURON': 'ALS inhibitor (Group 2)',
    'THIENCARBAZONE-METHYL': 'ALS inhibitor (Group 2)',
    'THIFENSULFURON': 'ALS inhibitor (Group 2)',
    'TRIASULFURON': 'ALS inhibitor (Group 2)',
    'TRIBENURON': 'ALS inhibitor (Group 2)',
    'TRIFLOXYSULFURON': 'ALS inhibitor (Group 2)',
    'TRIFLUSULFURON': 'ALS inhibitor (Group 2)',
    
    # EPSPS Inhibitors (Group 9) - 2 compounds
    'GLYPHOSATE': 'EPSPS Inhibitor (Group 9)',
    'SULFOSATE': 'EPSPS Inhibitor (Group 9)',
    
    # PPO inhibitors (Group 14) - 12 compounds
    'ACIFLUORFEN': 'PPO inhibitor (Group 14)',
    'CARFENTRAZONE-ETHYL': 'PPO inhibitor (Group 14)',
    'FLUMICLORAC': 'PPO inhibitor (Group 14)',
    'FLUMIOXAZIN': 'PPO inhibitor (Group 14)',
    'FLUTHIACET-METHYL': 'PPO inhibitor (Group 14)',
    'FOMESAFEN': 'PPO inhibitor (Group 14)',
    'LACTOFEN': 'PPO inhibitor (Group 14)',
    'OXADIAZON': 'PPO inhibitor (Group 14)',
    'OXYFLUORFEN': 'PPO inhibitor (Group 14)',
    'PYRAFLUFEN': 'PPO inhibitor (Group 14)',
    'SAFLUFENACIL': 'PPO inhibitor (Group 14)',
    'SULFENTRAZONE': 'PPO inhibitor (Group 14)',
    
    # Synthetic Auxins (Group 4) - 16 compounds
    '2,4-D': 'Synthetic Auxin (Group 4)',
    '2,4-DB': 'Synthetic Auxin (Group 4)',
    'AMINOPYRALID': 'Synthetic Auxin (Group 4)',
    'CHLORAMBEN': 'Synthetic Auxin (Group 4)',
    'CLOPYRALID': 'Synthetic Auxin (Group 4)',
    'DICAMBA': 'Synthetic Auxin (Group 4)',
    'DICHLORPROP': 'Synthetic Auxin (Group 4)',
    'FLORPYRAUXIFEN': 'Synthetic Auxin (Group 4)',
    'FLUROXYPYR': 'Synthetic Auxin (Group 4)',
    'HALAUXIFEN': 'Synthetic Auxin (Group 4)',
    'MCPA': 'Synthetic Auxin (Group 4)',
    'MCPB': 'Synthetic Auxin (Group 4)',
    'MECOPROP': 'Synthetic Auxin (Group 4)',
    'PICLORAM': 'Synthetic Auxin (Group 4)',
    'QUINCLORAC': 'Synthetic Auxin (Group 4)',
    'TRICLOPYR': 'Synthetic Auxin (Group 4)',
    
    # PS II Inhibitors - Group 5 (Ser264-binding site) - 15 compounds  
    'ATRAZINE': 'PS II Inhibitor (Group 5) - Triazine',
    'PROMETRYN': 'PS II Inhibitor (Group 5) - Triazine',
    'PROPAZINE': 'PS II Inhibitor (Group 5) - Triazine',
    'SIMAZINE': 'PS II Inhibitor (Group 5) - Triazine',
    'METRIBUZIN': 'PS II Inhibitor (Group 5) - Triazinone',
    'HEXAZINONE': 'PS II Inhibitor (Group 5) - Triazinone',
    'BROMACIL': 'PS II Inhibitor (Group 5) - Uracil',
    'TERBACIL': 'PS II Inhibitor (Group 5) - Uracil',
    'DIURON': 'PS II Inhibitor (Group 5) - Urea',
    'FLUOMETURON': 'PS II Inhibitor (Group 5) - Urea',
    'LINURON': 'PS II Inhibitor (Group 5) - Urea',
    'TEBUTHIURON': 'PS II Inhibitor (Group 5) - Urea',
    'DESMEDIPHAM': 'PS II Inhibitor (Group 5) - Phenylcarbamate',
    'PHENMEDIPHAM': 'PS II Inhibitor (Group 5) - Phenylcarbamate',
    'PROPANIL': 'PS II Inhibitor (Group 5) - Anilide',
    
    # PS II Inhibitors - Group 6 (His215-binding/alternative mechanisms) - 3 compounds
    'BENTAZONE': 'PS II Inhibitor (Group 6) - Benzothiadiazinone',
    'BROMOXYNIL': 'PS II Inhibitor (Group 6/24) - Nitrile (dual mode)',
    'DICHLOBENIL': 'PS II Inhibitor (Group 6) - Nitrile',
    
    # ACCase inhibitors (Group 1) - 8 compounds
    'CLETHODIM': 'ACCase inhibitor (Group 1) - DIM',
    'CLODINAFOP': 'ACCase inhibitor (Group 1) - FOP',
    'CYHALOFOP': 'ACCase inhibitor (Group 1) - FOP',
    'FLUAZIFOP': 'ACCase inhibitor (Group 1) - FOP',
    'PINOXADEN': 'ACCase inhibitor (Group 1) - DEN',
    'QUIZALOFOP': 'ACCase inhibitor (Group 1) - FOP',
    'SETHOXYDIM': 'ACCase inhibitor (Group 1) - DIM',
    'TRALKOXYDIM': 'ACCase inhibitor (Group 1) - DIM',
    
    # VLCFA Inhibitors (Group 15) - 12 compounds
    'ACETOCHLOR': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'ALACHLOR': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'DIMETHENAMID': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'DIMETHENAMID-P': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'EPTC': 'VLCFA Inhibitor (Group 15) - Thiocarbamate',
    'FLUFENACET': 'VLCFA Inhibitor (Group 15) - Oxyacetamide',
    'METOLACHLOR': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'METOLACHLOR-S': 'VLCFA Inhibitor (Group 15) - Chloroacetamide',
    'NAPROPAMIDE': 'VLCFA Inhibitor (Group 15) - Acetamide',
    'PYROXASULFONE': 'VLCFA Inhibitor (Group 15) - Isoxazoline',
    'THIOBENCARB': 'VLCFA Inhibitor (Group 15) - Thiocarbamate',
    'TRI-ALLATE': 'VLCFA Inhibitor (Group 15) - Thiocarbamate',
    
    # HPPD inhibitors (Group 27) - 6 compounds
    'BICYCLOPYRONE': 'HPPD inhibitor (Group 27)',
    'ISOXAFLUTOLE': 'HPPD inhibitor (Group 27)',
    'MESOTRIONE': 'HPPD inhibitor (Group 27)',
    'PYRASULFOTOLE': 'HPPD inhibitor (Group 27)',
    'TEMBOTRIONE': 'HPPD inhibitor (Group 27)',
    'TOPRAMEZONE': 'HPPD inhibitor (Group 27)',
    
    # Microtubule Assembly Inhibitors (Group 3) - 9 compounds
    'BENFLURALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'BUTRALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'DITHIOPYR': 'Microtubule Assembly Inhibitor (Group 3) - Pyridine',
    'ETHALFLURALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'ORYZALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'PENDIMETHALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'PRODIAMINE': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    'PROPYZAMIDE': 'Microtubule Assembly Inhibitor (Group 3) - Benzamide',
    'TRIFLURALIN': 'Microtubule Assembly Inhibitor (Group 3) - Dinitroaniline',
    
    # Other Herbicide MOAs - 15 compounds
    'CHLORPROPHAM': 'Mitosis Inhibitor (Group 18) - Carbamate',
    'CLOMAZONE': 'Carotenoid Biosynthesis Inhibitor (Group 13)',
    'DCPA': 'Mitosis Inhibitor (Group 3) - Benzoic acid',
    'DIFLUFENZOPYR': 'Auxin Transport Inhibitor (Group 19)',
    'DINOSEB': 'Uncoupler (Group 24) - Dinitrophenol',
    'DIQUAT': 'PS I Electron Diverter (Group 22) - Bipyridilium',
    'ETHOFUMESATE': 'Lipid synthesis inhibitor (Group 8)',
    'FLURIDONE': 'Carotenoid Biosynthesis Inhibitor (Group 12)',
    'GLUFOSINATE': 'Glutamine Synthetase Inhibitor (Group 10)',
    'INDAZIFLAM': 'Cellulose Biosynthesis Inhibitor (Group 29)',
    'ISOXABEN': 'Cellulose Biosynthesis Inhibitor (Group 21)',
    'MSMA': 'Unknown MOA (Group 17) - Organoarsenic',
    'NORFLURAZON': 'Carotenoid Biosynthesis Inhibitor (Group 12)',
    'PARAQUAT': 'PS I Electron Diverter (Group 22) - Bipyridilium',
    'PELARGONIC': 'Cell membrane disruptor (Group 26)'
}# -------------------------------------------------------------------------

# Aliases for common variants
ALIASES = {
    "S-METOLACHLOR": "METOLACHLOR-S",
    "TRIALLATE": "TRI-ALLATE",
    "GLUFOSINATE-AMMONIUM": "GLUFOSINATE",
    "PYRAFLUFEN-ETHYL": "PYRAFLUFEN",
    "HALAUXIFEN-METHYL": "HALAUXIFEN",
    "FLORPYRAUXIFEN-BENZYL": "FLORPYRAUXIFEN",
    "BROMOXYNIL OCTANOATE": "BROMOXYNIL",
    "DICHLORPROP-P": "DICHLORPROP",
}
def log(msg: str):
    if VERBOSE:
        print(msg, file=sys.stderr)

def normalize_compound(name: str) -> str:
    if pd.isna(name):
        return name
    s = str(name).upper().strip()
    base = s.split(" ")[0]  # drop salt/formulation suffixes
    return ALIASES.get(s, ALIASES.get(base, s))

def pick_value(kind: str, df: pd.DataFrame) -> pd.Series:
    if kind == "high":
        return pd.to_numeric(df["EPEST_HIGH_KG"], errors="coerce")
    if kind == "low":
        return pd.to_numeric(df["EPEST_LOW_KG"], errors="coerce")
    if kind == "mean":
        return (pd.to_numeric(df["EPEST_LOW_KG"], errors="coerce") +
                pd.to_numeric(df["EPEST_HIGH_KG"], errors="coerce")) / 2.0
    raise ValueError("VALUE_KIND must be 'high' | 'low' | 'mean'")

# =========================
# Sites → county FIPS
# =========================
def load_sites_with_fips(sites_csv: str, county_shp: str) -> pd.DataFrame:
    counties = gpd.read_file(county_shp)
    if "GEOID" in counties.columns:
        counties = counties.rename(columns={"GEOID": "fips"})
    assert "fips" in counties.columns, "County shapefile must include GEOID or fips"
    counties = counties[["fips", "geometry"]].to_crs(4326)

    sites = pd.read_csv(sites_csv)
    for c in ["decimalLongitude", "decimalLatitude", "year"]:
        if c not in sites.columns:
            raise ValueError(f"Sites CSV missing required column: {c}")

    pts = gpd.GeoDataFrame(
        sites.copy(),
        geometry=[Point(xy) for xy in zip(sites["decimalLongitude"], sites["decimalLatitude"])],
        crs="EPSG:4326"
    )

    joined = gpd.sjoin(pts, counties, how="left", predicate="within")
    out = pd.DataFrame(joined.drop(columns=["geometry", "index_right"]))
    out["fips"] = out["fips"].astype(str).str.zfill(5)
    out["year"] = pd.to_numeric(out["year"], errors="coerce").astype("Int64")
    return out

# =========================
# EPest loading (URL or local)
# =========================
def _source_for_year(y: int) -> tuple[str, str]:
    """
    Return ("local", path) if overridden; else ("url", url)
    """
    y = int(y)
    if y in LOCAL_BY_YEAR:
        return "local", LOCAL_BY_YEAR[y]
    return "url", EPEST_URL.format(year=y)

def load_epest_year(year: int, value_kind: str) -> pd.DataFrame:
    y = int(year)
    src_kind, src = _source_for_year(y)
    log(f"[EPest] Loading {y} from {src_kind.upper()}: {src}")

    usecols = [
        "COMPOUND", "YEAR", "STATE_FIPS_CODE", "COUNTY_FIPS_CODE",
        "EPEST_LOW_KG", "EPEST_HIGH_KG"
    ]

    try:
        # pd.read_csv supports both URLs and local paths
        df = pd.read_csv(src, sep="\t", usecols=usecols, dtype=str, low_memory=False)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Local EPest file for {y} not found at: {src}\n"
            f"Update LOCAL_BY_YEAR[{y}] to a valid path."
        ) from e

    # Standardize FIPS and YEAR
    df["STATE_FIPS_CODE"]  = df["STATE_FIPS_CODE"].str.zfill(2)
    df["COUNTY_FIPS_CODE"] = df["COUNTY_FIPS_CODE"].str.zfill(3)
    df["county_fips"] = df["STATE_FIPS_CODE"] + df["COUNTY_FIPS_CODE"]
    df["YEAR"] = pd.to_numeric(df["YEAR"], errors="coerce").astype("Int64")

    # Map compound -> MOA
    df["COMPOUND_NORM"] = df["COMPOUND"].map(normalize_compound)
    df["MOA"] = df["COMPOUND_NORM"].map(epest_herbicide_moa).fillna("Other/Unmapped")

    # Pick value
    df["kg"] = pick_value(value_kind, df)

    # Aggregate county × year × MOA
    agg = (df.groupby(["YEAR", "county_fips", "MOA"], dropna=False)["kg"]
             .sum(min_count=1).reset_index())
    return agg

def aggregate_epest_to_wide(years_needed: list[int], value_kind: str) -> pd.DataFrame:
    frames = []
    for y in years_needed:
        frames.append(load_epest_year(y, value_kind))
    tall = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["YEAR", "county_fips", "MOA", "kg"]
    )
    wide = (tall.pivot_table(index=["YEAR", "county_fips"], columns="MOA", values="kg", aggfunc="sum")
                 .fillna(0.0)
                 .reset_index())
    wide.columns.name = None
    return wide

# =========================
# Build final table
# =========================
def build_sites_with_action_sites() -> pd.DataFrame:
    sites = load_sites_with_fips(SITES_CSV, COUNTY_SHP)

    # Important: cast to int so override checks (e.g., 2018) work
    years_needed = sorted(int(y) for y in sites["year"].dropna().unique())
    log(f"[Sites] Years needed: {years_needed}")

    epest_wide = aggregate_epest_to_wide(years_needed, VALUE_KIND)

    merged = sites.merge(
        epest_wide,
        left_on=["year", "fips"],
        right_on=["YEAR", "county_fips"],
        how="left"
    ).drop(columns=["YEAR", "county_fips"])

    moa_cols = [c for c in merged.columns if c not in sites.columns]
    merged[moa_cols] = merged[moa_cols].fillna(0.0)

    if OUT_PATH.lower().endswith(".parquet"):
        merged.to_parquet(OUT_PATH, index=False)
    else:
        merged.to_csv(OUT_PATH, index=False)

    return merged

# =========================
# MAIN
# =========================
if __name__ == "__main__":
    try:
        df_out = build_sites_with_action_sites()
    except Exception as e:
        # Surface a clear, single-line message, then re-raise for full traceback if needed
        print(f"\nERROR: {e}\n", file=sys.stderr)
        raise

    print(f"Done. Wrote: {OUT_PATH}")
    base_cols = {"decimalLongitude", "decimalLatitude", "year", "fips"}
    new_cols = [c for c in df_out.columns if c not in base_cols]
    print("New columns (sample):", new_cols[:8])