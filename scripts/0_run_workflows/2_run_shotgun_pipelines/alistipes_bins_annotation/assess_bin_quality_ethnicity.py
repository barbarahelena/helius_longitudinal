#!/usr/bin/env python3
"""
Assess bin quality with stricter thresholds (completeness >= 80%, contamination < 5%).
For passing bins, traces back to the source sample and merges with clinical data
to report the ethnic group distribution per taxon.

Source sample is determined as the sample with maximum depth coverage in
the bin_summary.tsv files (depth columns: "Depth {sample_id}").
"""

import re
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
BIN_SUMMARY_FILES = [
    BASE_DIR / "mag_results" / "results_batch1" / "GenomeBinning" / "bin_summary.tsv",
    BASE_DIR / "mag_results" / "results_batch2" / "GenomeBinning" / "bin_summary.tsv",
    BASE_DIR / "mag_results" / "results_batch3" / "GenomeBinning" / "bin_summary.tsv",
]

# Clinical data: semicolon-delimited CSV with columns Sample_ID, Heliusnr, Ethnicity
CLINICAL_CSV = Path(__file__).parents[5] / "data" / "clinicaldata" / "clinical_data_long_03_07_2023_N3443.csv"

COMPLETENESS_THRESHOLD = 80.0
CONTAMINATION_THRESHOLD = 5.0


def extract_heliusnr(sample_name: str) -> str | None:
    """Extract numeric subject ID from sample names like HELIBA_100245 or HELIFU_100245."""
    m = re.search(r"HELI\w+[_.](\d+)", sample_name)
    return m.group(1) if m else None


# --- Load CheckM2 summaries ---
checkm2 = pd.concat(
    [pd.read_csv(f, sep="\t") for f in CHECKM2_FILES if f.exists()],
    ignore_index=True,
)
checkm2["Name"] = checkm2["Name"].str.strip()
print(f"Total bins in CheckM2 summaries: {len(checkm2)}\n")

# --- Load bin_summary.tsv files and build bin -> source_sample mapping ---
# Source sample = sample with highest depth coverage for each bin
depth_dfs = []
for f in BIN_SUMMARY_FILES:
    if not f.exists():
        print(f"WARNING: bin_summary.tsv not found: {f}")
        continue
    df = pd.read_csv(f, sep="\t", low_memory=False)
    depth_cols = [c for c in df.columns if c.startswith("Depth ")]
    if not depth_cols or "bin" not in df.columns:
        print(f"WARNING: expected columns not found in {f}")
        continue
    bin_col = df["bin"].str.replace(r"\.fa$", "", regex=True).str.strip()
    depth_matrix = df[depth_cols].copy()
    depth_matrix.index = bin_col
    depth_dfs.append(depth_matrix)

bin_to_sample: dict[str, str] = {}
if depth_dfs:
    all_depths = pd.concat(depth_dfs)
    # Remove duplicate bin entries (keep first)
    all_depths = all_depths[~all_depths.index.duplicated(keep="first")]
    for bin_id, row in all_depths.iterrows():
        if row.max() > 0:
            best_col = row.idxmax()  # e.g., "Depth HELIBA_100245"
            bin_to_sample[str(bin_id)] = best_col.replace("Depth ", "").strip()
else:
    print("WARNING: No bin_summary.tsv files loaded — ethnicity linkage will be skipped.\n")

# --- Load clinical data ---
clinical: pd.DataFrame | None = None
if CLINICAL_CSV.exists():
    clinical = pd.read_csv(CLINICAL_CSV, sep=";")
    # Keep one row per Heliusnr (use baseline row if duplicated)
    clinical["Heliusnr"] = clinical["Heliusnr"].astype(str)
    clinical = clinical.drop_duplicates(subset="Heliusnr", keep="first")
    print(f"Clinical data loaded: {len(clinical)} subjects\n")
else:
    print(f"WARNING: Clinical CSV not found at {CLINICAL_CSV} — ethnicity linkage will be skipped.\n")

# --- Assess each taxon group ---
for taxon, samplesheet_path in SAMPLESHEETS.items():
    if not samplesheet_path.exists():
        print(f"WARNING: Samplesheet not found: {samplesheet_path}\n")
        continue

    ss = pd.read_csv(samplesheet_path)
    id_col = "sample" if "sample" in ss.columns else "bin"
    bin_names = set(ss[id_col].str.strip())

    taxon_qc = checkm2[checkm2["Name"].isin(bin_names)].copy()
    total = len(taxon_qc)
    missing = len(bin_names) - total

    hq = taxon_qc[
        (taxon_qc["Completeness"] >= COMPLETENESS_THRESHOLD) &
        (taxon_qc["Contamination"] < CONTAMINATION_THRESHOLD)
    ].copy()
    n_hq = len(hq)
    pct_hq = (n_hq / total * 100) if total > 0 else 0.0

    print(f"=== {taxon} ===")
    print(f"  Bins in samplesheet:     {len(bin_names)}")
    print(f"  Found in CheckM2:        {total}")
    if missing:
        not_found = bin_names - set(checkm2["Name"])
        print(f"  WARNING — not in CheckM2: {sorted(not_found)}")
    print(f"  Passing QC (>={COMPLETENESS_THRESHOLD}% complete, <{CONTAMINATION_THRESHOLD}% contam): "
          f"{n_hq} / {total} ({pct_hq:.1f}%)")

    # Distribution stats for non-passing bins
    lq = taxon_qc[~taxon_qc.index.isin(hq.index)]
    if len(lq) > 0:
        fail_c  = lq[(lq["Completeness"] < COMPLETENESS_THRESHOLD) & (lq["Contamination"] < CONTAMINATION_THRESHOLD)]
        fail_co = lq[(lq["Completeness"] >= COMPLETENESS_THRESHOLD) & (lq["Contamination"] >= CONTAMINATION_THRESHOLD)]
        fail_b  = lq[(lq["Completeness"] < COMPLETENESS_THRESHOLD) & (lq["Contamination"] >= CONTAMINATION_THRESHOLD)]
        print(f"  Non-passing ({len(lq)}): "
              f"low completeness only={len(fail_c)}, "
              f"high contamination only={len(fail_co)}, "
              f"both={len(fail_b)}")

    # --- Ethnicity breakdown ---
    if bin_to_sample and clinical is not None and n_hq > 0:
        hq_bins = hq["Name"].tolist()
        source_samples = [bin_to_sample.get(b) for b in hq_bins]
        heliusnrs = [extract_heliusnr(s) for s in source_samples if s is not None]
        heliusnrs_found = [h for h in heliusnrs if h is not None]

        n_no_sample = sum(1 for s in source_samples if s is None)
        n_no_heliusnr = len([s for s in source_samples if s is not None]) - len(heliusnrs_found)

        if n_no_sample:
            print(f"  WARNING: {n_no_sample} HQ bin(s) not found in depth data")
        if n_no_heliusnr:
            print(f"  WARNING: {n_no_heliusnr} sample(s) could not be parsed for Heliusnr")

        merged = (
            pd.DataFrame({"Heliusnr": heliusnrs_found})
            .merge(clinical[["Heliusnr", "Ethnicity"]], on="Heliusnr", how="left")
        )
        eth_counts = merged["Ethnicity"].value_counts(dropna=False)
        print(f"\n  Ethnicity of subjects with HQ bins (n unique subjects = "
              f"{merged['Heliusnr'].nunique()}, total HQ bins = {n_hq}):")
        for eth, count in eth_counts.items():
            print(f"    {eth}: {count} bins")
    print()

# --- Save detailed HQ bin table ---
all_hq_rows = []
for taxon, samplesheet_path in SAMPLESHEETS.items():
    if not samplesheet_path.exists():
        continue
    ss = pd.read_csv(samplesheet_path)
    id_col = "sample" if "sample" in ss.columns else "bin"
    bin_names = set(ss[id_col].str.strip())
    taxon_qc = checkm2[checkm2["Name"].isin(bin_names)].copy()
    hq = taxon_qc[
        (taxon_qc["Completeness"] >= COMPLETENESS_THRESHOLD) &
        (taxon_qc["Contamination"] < CONTAMINATION_THRESHOLD)
    ].copy()
    hq["Taxon"] = taxon
    hq["Source_sample"] = hq["Name"].map(bin_to_sample)
    hq["Heliusnr"] = hq["Source_sample"].apply(
        lambda s: extract_heliusnr(s) if pd.notna(s) else None
    )
    all_hq_rows.append(hq)

if all_hq_rows:
    hq_all = pd.concat(all_hq_rows, ignore_index=True)
    if clinical is not None:
        hq_all = hq_all.merge(
            clinical[["Heliusnr", "Ethnicity", "Sex"]],
            on="Heliusnr", how="left"
        )
    out_path = BASE_DIR / "data" / "bin_quality_hq_ethnicity.tsv"
    cols = ["Taxon", "Name", "Completeness", "Contamination", "Source_sample",
            "Heliusnr", "Ethnicity", "Sex"]
    cols_present = [c for c in cols if c in hq_all.columns]
    hq_all[cols_present].to_csv(out_path, sep="\t", index=False)
    print(f"Detailed HQ bin table saved to: {out_path}")
