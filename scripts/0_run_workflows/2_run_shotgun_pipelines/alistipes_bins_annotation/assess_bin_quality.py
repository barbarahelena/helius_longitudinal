#!/usr/bin/env python3
"""
Assess bin quality for Alistipes, Parabacteroides, and Senegalimassilia samplesheets.
Filters for bins with completeness > 70% and contamination < 10%.
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

# --- Load and merge all CheckM2 summaries ---
checkm2_dfs = []
for f in CHECKM2_FILES:
    if f.exists():
        df = pd.read_csv(f, sep="\t")
        checkm2_dfs.append(df)
    else:
        print(f"WARNING: CheckM2 file not found: {f}")

checkm2 = pd.concat(checkm2_dfs, ignore_index=True)
print(f"Total bins in CheckM2 summaries: {len(checkm2)}\n")

# --- Assess each taxon group ---
results = []

for taxon, samplesheet_path in SAMPLESHEETS.items():
    ss = pd.read_csv(samplesheet_path)
    # Support both 'sample' and 'bin' column names
    id_col = "sample" if "sample" in ss.columns else "bin"
    bin_names = set(ss[id_col].str.strip())

    # Match bins in checkm2 (Name column, strip whitespace)
    taxon_qc = checkm2[checkm2["Name"].str.strip().isin(bin_names)].copy()

    total = len(taxon_qc)
    missing = len(bin_names) - total
    hq = taxon_qc[
        (taxon_qc["Completeness"] > COMPLETENESS_THRESHOLD) &
        (taxon_qc["Contamination"] < CONTAMINATION_THRESHOLD)
    ]
    n_hq = len(hq)
    pct_hq = (n_hq / total * 100) if total > 0 else 0.0

    results.append({
        "Taxon":            taxon,
        "Total_bins":       len(bin_names),
        "Found_in_CheckM2": total,
        "Missing_from_CheckM2": missing,
        f"HQ (>={COMPLETENESS_THRESHOLD}% complete, <{CONTAMINATION_THRESHOLD}% contam)": n_hq,
        "% HQ":             round(pct_hq, 1),
    })

    # Detailed summary
    print(f"=== {taxon} ===")
    print(f"  Bins in samplesheet:     {len(bin_names)}")
    print(f"  Found in CheckM2:        {total}")
    if missing:
        not_found = bin_names - set(checkm2["Name"].str.strip())
        print(f"  WARNING - not in CheckM2: {sorted(not_found)}")
    print(f"  HQ bins (>{COMPLETENESS_THRESHOLD}% complete, <{CONTAMINATION_THRESHOLD}% contam): {n_hq} / {total} ({pct_hq:.1f}%)")

    # Distribution of completeness for non-HQ bins
    lq = taxon_qc[~taxon_qc.index.isin(hq.index)]
    if len(lq) > 0:
        print(f"  --- Non-HQ bins ({len(lq)}) ---")
        print(f"    Completeness: mean={lq['Completeness'].mean():.1f}%, "
              f"min={lq['Completeness'].min():.1f}%, max={lq['Completeness'].max():.1f}%")
        print(f"    Contamination: mean={lq['Contamination'].mean():.1f}%, "
              f"min={lq['Contamination'].min():.1f}%, max={lq['Contamination'].max():.1f}%")
        # Show bins failing each criterion
        fail_complete = lq[lq["Completeness"] <= COMPLETENESS_THRESHOLD]
        fail_contam   = lq[lq["Contamination"] >= CONTAMINATION_THRESHOLD]
        print(f"    Failing completeness only: {len(fail_complete[fail_complete['Contamination'] < CONTAMINATION_THRESHOLD])}")
        print(f"    Failing contamination only: {len(fail_contam[fail_contam['Completeness'] > COMPLETENESS_THRESHOLD])}")
        print(f"    Failing both: {len(lq[(lq['Completeness'] <= COMPLETENESS_THRESHOLD) & (lq['Contamination'] >= CONTAMINATION_THRESHOLD)])}")
    print()

# --- Summary table ---
summary_df = pd.DataFrame(results)
print("=== SUMMARY TABLE ===")
print(summary_df.to_string(index=False))

# Save summary
out_path = BASE_DIR / "data" / "bin_quality_summary.tsv"
summary_df.to_csv(out_path, sep="\t", index=False)
print(f"\nSummary saved to: {out_path}")
