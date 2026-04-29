#!/usr/bin/env python3
"""
Filter bin_summary.tsv for bins classified as Parabacteroides distasonis
and save results to a new CSV file in the same directory.
"""

import pandas as pd
from pathlib import Path

input_file = Path(__file__).parent.parent / "mag_results/results_batch3/GenomeBinning/bin_summary.tsv"
output_file = Path(__file__).parent.parent / "mag_results/results_batch3/GenomeBinning/parabacteroides_distasonis_bins.csv"

df = pd.read_csv(input_file, sep="\t", low_memory=False)

# Search for Parabacteroides distasonis in all taxonomy/classification columns
taxonomy_cols = [col for col in df.columns if any(
    kw in col.lower() for kw in ["classification", "taxonomy"]
)]

mask = df[taxonomy_cols].apply(
    lambda col: col.astype(str).str.contains("Parabacteroides distasonis", case=False, na=False)
).any(axis=1)

filtered = df[mask]

print(f"Found {len(filtered)} bins classified as Parabacteroides distasonis")
print(f"Taxonomy columns searched: {taxonomy_cols}")

if len(filtered) > 0:
    print("\nBin names:")
    for name in filtered["Name"]:
        print(f"  {name}")

filtered.to_csv(output_file, index=False)
print(f"\nSaved to: {output_file}")
