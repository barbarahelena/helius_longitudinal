#!/bin/bash
# Count CDS features per Bakta TSV and write a summary table.
# Usage: bash count_cds.sh [input_dir] [output_file]

indir="${1:-.}"
out="${2:-cds_summary.tsv}"

printf "ID\tn_cds\n" > "$out"

for f in "$indir"/*.tsv; do
    [[ "$(basename "$f")" == "$(basename "$out")" ]] && continue
    id=$(basename "$f" .tsv)
    n=$(awk -F'\t' '!/^#/ && $2 == "cds"' "$f" | wc -l)
    printf "%s\t%s\n" "$id" "$n" >> "$out"
done

echo "Wrote $(($(wc -l < "$out") - 1)) entries to $out"
