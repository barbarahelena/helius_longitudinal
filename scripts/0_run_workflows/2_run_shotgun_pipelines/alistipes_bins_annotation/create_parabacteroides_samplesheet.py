#!/usr/bin/env python3
"""
Creates a samplesheet for all Parabacteroides distasonis bins across 3 assembly batches.
Output: data/parabacteroides_samplesheet.csv with columns: bin, path
"""

import csv
import os

BASE = "mag_results"
BATCHES = ["results_batch1", "results_batch2", "results_batch3"]
OUTPUT = "data/parabacteroides_samplesheet.csv"

rows = []
missing = []

for batch in BATCHES:
    csv_path = os.path.join(BASE, batch, "GenomeBinning", "parabacteroides_distasonis_bins.csv")
    bins_dir = os.path.join(BASE, batch, "GenomeBinning", "DASTool", "bins")

    if not os.path.exists(csv_path):
        print(f"WARNING: CSV not found: {csv_path}")
        continue

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bin_filename = row["bin"].strip()
            # Strip .fa suffix if present to get the bin name
            bin_name = bin_filename.replace(".fa", "")
            bin_path = os.path.join(bins_dir, bin_filename)

            if os.path.exists(bin_path):
                rows.append({"bin": bin_name, "path": bin_path})
            else:
                missing.append((batch, bin_filename, bin_path))

with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["bin", "path"])
    writer.writeheader()
    writer.writerows(rows)

print(f"Written {len(rows)} bins to {OUTPUT}")

if missing:
    print(f"\nWARNING: {len(missing)} bins not found in DASTool/bins:")
    for batch, name, path in missing:
        print(f"  [{batch}] {name} -> {path}")
