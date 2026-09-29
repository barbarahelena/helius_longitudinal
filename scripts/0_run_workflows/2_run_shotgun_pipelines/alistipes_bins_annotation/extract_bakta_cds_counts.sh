#!/bin/bash
# Per-bin CDS counts from the Bakta annotation, for use as the denominator in
# the VFDB and KEGG gene-content proportions.
#
# The count is taken from each bin's Bakta protein FASTA, so it is exactly the
# set of proteins that DIAMOND searched against VFDB — numerator and
# denominator then come from the same gene prediction.
#
# Run on the login node, then copy the output next to the other annotation data:
#   bash extract_bakta_cds_counts.sh
#   # -> $BASE_DIR/bakta_cds_counts.tsv
#   # copy to data/shotgun/alistipes_annotation/bakta_cds_counts.tsv in the repo
#
# Locus prefixes are read from the FASTA headers (e.g. >FNBCDG_00001 -> FNBCDG),
# not from filenames, so the layout of the Bakta output directory does not matter.

set -euo pipefail

BASE_DIR="${BASE_DIR:-/projects/prjs1784/heliuspaired}"
BAKTA_DIR="${BAKTA_DIR:-${BASE_DIR}/annotation_alistipes/bakta}"
OUT="${OUT:-${BASE_DIR}/bakta_cds_counts.tsv}"

[[ -d "$BAKTA_DIR" ]] || {
  echo "ERROR: Bakta output directory not found: ${BAKTA_DIR}" >&2
  echo "       Set BAKTA_DIR=/path/to/bakta (the annotation pipeline's outdir)." >&2
  exit 1
}

mapfile -t faas < <(find "$BAKTA_DIR" -name '*.faa' | sort)
[[ ${#faas[@]} -gt 0 ]] || {
  echo "ERROR: no .faa files under ${BAKTA_DIR}" >&2
  exit 1
}
echo "Bakta protein FASTAs found: ${#faas[@]}"

printf 'locus_prefix\tcds_count\tsource_file\n' > "$OUT"
for f in "${faas[@]}"; do
  n=$(grep -c '^>' "$f" || true)
  prefix=$(head -1 "$f" | sed -E 's/^>([A-Za-z0-9]+)_[0-9]+.*/\1/')
  if [[ -z "$prefix" || "$prefix" == ">"* ]]; then
    echo "WARNING: could not read a locus prefix from ${f}, skipping" >&2
    continue
  fi
  printf '%s\t%s\t%s\n' "$prefix" "$n" "$(basename "$f")" >> "$OUT"
done

n_rows=$(( $(wc -l < "$OUT") - 1 ))
n_uniq=$(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l)
echo "Wrote ${n_rows} rows (${n_uniq} distinct locus prefixes) to ${OUT}"
if [[ "$n_rows" -ne "$n_uniq" ]]; then
  echo "WARNING: some locus prefixes appear more than once — check for duplicate Bakta runs" >&2
fi
echo
echo "CDS count summary:"
tail -n +2 "$OUT" | cut -f2 | sort -n | awk '{a[NR]=$1} END {printf "  min %d  median %d  max %d\n", a[1], a[int((NR+1)/2)], a[NR]}'
