#!/usr/bin/env python3
"""
Filter samplesheets to retain only bins passing QC thresholds:
  - Completeness > 70%
  - Contamination < 10%

Outputs filtered samplesheets alongside the originals with a _filtered suffix.
"""

import pandas as pd
from pathlib import Path

# --- Config ---
BASE_DIR = Path("/projects/prjs1784/heliuspaired")
SAMPLESHEETS = {
    "Alistipes":        BASE_DIR / "data" / "alistipes_samplesheet.csv",
    "Parabacteroides":  BASE_DIR / "data" / "parabacteroides_samplesheet.csv",
    "Senegalimassilia": BASE_DIR / "data" / "senegalimassilia_samplesheet.csv",
}
CHECKM2_FILES = [
    BASE_DIR / "mag_results" / "summaries" / "batch1" / "checkm2_summary.tsv",
    BASE_DIR / "mag_results" / "summaries" / "batch2" / "checkm2_summary.tsv",
    BASE_DIR / "mag_results" / "summaries" / "batch3" / "checkm2_summary.tsv",
]
COMPLETENESS_THRESHOLD = 70.0
CONTAMINATION_THRESHOLD = 10.0

# --- Load CheckM2 ---
checkm2 = pd.concat(
    [pd.read_csv(f, sep="\t") for f in CHECKM2_FILES if f.exists()],
    ignore_index=True
)
checkm2["Name"] = checkm2["Name"].str.strip()

hq_bins = set(
    checkm2.loc[
        (checkm2["Completeness"] > COMPLETENESS_THRESHOLD) &
        (checkm2["Contamination"] < CONTAMINATION_THRESHOLD),
        "Name"
    ]
)

# --- Filter each samplesheet ---
for taxon, ss_path in SAMPLESHEETS.items():
    ss = pd.read_csv(ss_path)
    id_col = "sample" if "sample" in ss.columns else "bin"

    before = len(ss)
    ss_filtered = ss[ss[id_col].str.strip().isin(hq_bins)].copy()
    after = len(ss_filtered)
    removed = before - after

    out_path = ss_path.with_name(ss_path.stem + "_filtered.csv")
    ss_filtered.to_csv(out_path, index=False)

    print(f"{taxon}: {before} -> {after} bins ({removed} removed)")
    if removed > 0:
        removed_bins = ss[~ss[id_col].str.strip().isin(hq_bins)][id_col].tolist()
        for b in removed_bins:
            row = checkm2[checkm2["Name"] == b]
            if not row.empty:
                c = row.iloc[0]["Completeness"]
                cont = row.iloc[0]["Contamination"]
                print(f"  REMOVED: {b}  (completeness={c:.1f}%, contamination={cont:.1f}%)")
            else:
                print(f"  REMOVED: {b}  (not found in CheckM2)")
    print(f"  -> Written to {out_path}\n")
