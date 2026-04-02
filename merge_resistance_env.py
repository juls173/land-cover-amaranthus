
import pandas as pd
from pathlib import Path

# --- file paths ---
ENV_PATH = "sites_with_action_sites.csv"
GEN_PATH = "genetic_data\\resistancefreqs_herbariumandcontemp.txt"     # tab-separated
OUT_PATH = "resistance_gf_input.csv"

gen_raw = pd.read_csv(GEN_PATH, sep="\t")


env_raw = pd.read_csv(ENV_PATH)                 # if TSV: pd.read_csv(ENV_PATH, sep="\t")
gen_raw = pd.read_csv(GEN_PATH, sep="\t")

# --- unify keys: env uses SeqNames, genetics uses Sample ---
gen_raw = gen_raw.rename(columns={"Sample": "SeqNames"})

env_index = env_raw.columns.get_loc("year")
env_cols = ["decimalLongitude", "decimalLatitude"] + list(env_raw.columns[env_index:])

# 1) Verify ONLY those env columns are identical within SeqNames
cols_here = [c for c in env_cols if c in env_raw.columns]
varying = []
for name, g in env_raw.groupby("SeqNames"):
    if not g[cols_here].apply(lambda col: col.eq(col.iloc[0]).all()).all():
        varying.append(name)
if varying:
    raise ValueError(f"Environmental columns vary within these SeqNames (first 10): {varying[:10]}")

# 2) Deduplicate env to one row per SeqNames (keep env cols; coords optional)
env_dedup_env = env_raw.groupby("SeqNames", as_index=False)[cols_here].first()

coord_candidates = ["lon","lat","longitude","latitude","x","y"]
coord_cols = [c for c in coord_candidates if c in env_raw.columns]
if coord_cols:
    env_coords = env_raw.groupby("SeqNames", as_index=False)[coord_cols].mean(numeric_only=True)
    env_dedup = env_dedup_env.merge(env_coords, on="SeqNames", how="left")
else:
    env_dedup = env_dedup_env

snp_cols = [c for c in gen_raw.columns if c != "SeqNames"]
if not snp_cols:
    raise ValueError("No SNP columns.")

gen_raw[snp_cols] = gen_raw[snp_cols].apply(pd.to_numeric, errors="coerce")
gen_by_pop = gen_raw.groupby("SeqNames", as_index=False)[snp_cols].mean(numeric_only=True)

print(f"\nSample SeqNames from gen_by_pop: {list(gen_by_pop['SeqNames'].head(10))}")
print(f"Sample SeqNames from env_dedup: {list(env_dedup['SeqNames'].head(10))}")
print(f"\nOverlapping SeqNames: {set(gen_by_pop['SeqNames']) & set(env_dedup['SeqNames'])}")
# 4) Merge (no SNP dropping)
merged = env_dedup.merge(gen_by_pop, on="SeqNames", how="inner")
merged.to_csv(OUT_PATH, index=False)

print(f"Merged shape: {merged.shape}")
print(f"Kept ALL DroughtAnc_ SNPs: {len(snp_cols)}")
print(f"Wrote: {Path(OUT_PATH).resolve()}")