#!/usr/bin/env python3
"""
Build a translation sheet mapping bin_name -> locus_tag_prefix
from bakta GFF3 annotation output.

Usage:
    python make_translation_sheet.py
"""

import re
from pathlib import Path

BAKTA_DIRS = {
    "parabacteroides":  Path("/projects/prjs1784/heliuspaired/annotation_parabacteroides/bakta"),
    "senegalimassilia": Path("/projects/prjs1784/heliuspaired/annotation_senegalimassilia/bakta"),
}

for taxon, bakta_dir in BAKTA_DIRS.items():
    if not bakta_dir.exists():
        print(f"Skipping {taxon}: {bakta_dir} does not exist yet")
        continue

    gff3_files = sorted(bakta_dir.glob("*.gff3"))
    if not gff3_files:
        print(f"Skipping {taxon}: no GFF3 files found in {bakta_dir}")
        continue

    rows = []
    for gff3 in gff3_files:
        bin_name = gff3.stem  # filename without extension
        prefix = None
        with open(gff3) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                m = re.search(r'locus_tag=([^;_]+)', line)
                if m:
                    # prefix is everything before the first underscore+digits
                    tag = m.group(1)
                    pm = re.match(r'^([A-Z]+)', tag)
                    prefix = pm.group(1) if pm else tag
                    break
        if prefix is None:
            print(f"  WARNING: no locus_tag found for {bin_name}")
            prefix = "NA"
        rows.append((bin_name, prefix))

    out_path = Path("/projects/prjs1784/heliuspaired/data") / f"{taxon}_translation_sheet.tsv"
    with open(out_path, "w") as out:
        out.write("bin_name\tlocus_tag_prefix\n")
        for bin_name, prefix in rows:
            out.write(f"{bin_name}\t{prefix}\n")

    print(f"Written {len(rows)} entries to {out_path}")
