#!/usr/bin/env python3
"""
Build the phylogenomics samplesheet (sample,gff3) for Alistipes bins with
CheckM2 completeness >= 70% and contamination < 10%.

Input:  helius_longitudinal/data/alistipes_samplesheet.csv  (bin,path or sample,fasta)
Output: data/alistipes_bakta_samplesheet_c70.csv            (sample,gff3)
Bins without an existing Bakta GFF3 are reported and dropped.
"""
import pandas as pd
from pathlib import Path

BASE_DIR = Path("/projects/prjs1784/heliuspaired")
IN_SHEET = BASE_DIR / "helius_longitudinal" / "data" / "alistipes_samplesheet.csv"
OUT_SHEET = BASE_DIR / "data" / "alistipes_bakta_samplesheet_c70.csv"
CHECKM2_FILES = [BASE_DIR / "mag_results" / "summaries" / f"batch{i}" / "checkm2_summary.tsv" for i in (1, 2, 3)]
GFF_DIR = "annotation_alistipes/bakta"   # relative to BASE_DIR, as the pipeline expects
COMPLETENESS_MIN = 70.0
CONTAMINATION_MAX = 10.0

checkm2 = pd.concat([pd.read_csv(f, sep="\t") for f in CHECKM2_FILES], ignore_index=True)
checkm2["Name"] = checkm2["Name"].str.strip()

bins = pd.read_csv(IN_SHEET)
# accept either header: bin,path or sample,fasta
id_col = "bin" if "bin" in bins.columns else "sample"
bins["bin"] = bins[id_col].str.strip()
merged = bins.merge(checkm2[["Name", "Completeness", "Contamination"]],
                    left_on="bin", right_on="Name", how="left")

missing = merged[merged["Completeness"].isna()]
if len(missing):
    raise SystemExit(f"{len(missing)} bins not found in CheckM2: {missing['bin'].tolist()}")

keep = merged[(merged["Completeness"] >= COMPLETENESS_MIN) & (merged["Contamination"] < CONTAMINATION_MAX)].copy()
keep["gff3"] = keep["bin"].map(lambda b: f"{GFF_DIR}/{b}.gff3")
no_gff = keep[~keep["gff3"].map(lambda p: (BASE_DIR / p).exists())]
if len(no_gff):
    print(f"WARNING: dropping {len(no_gff)} bins without a Bakta GFF3: {no_gff['bin'].tolist()}")
    keep = keep.drop(no_gff.index)

keep[["bin", "gff3"]].rename(columns={"bin": "sample"}).to_csv(OUT_SHEET, index=False)
print(f"{len(bins)} bins -> {len(keep)} with completeness >= {COMPLETENESS_MIN}% and contamination < {CONTAMINATION_MAX}%")
print(f"Written to {OUT_SHEET}")
