#!/bin/bash
#SBATCH --job-name=instrain_ap_between
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --array=1-130%20
#SBATCH --output=logs/instrain_between_%A_%a.out
#SBATCH --error=logs/instrain_between_%A_%a.err

# Alistipes putredinis strain retention: BETWEEN-person background comparison.
#
# Why: 1_instrain_strain_retention.R (in the alistipes chapter) compares each
# participant's own baseline vs follow-up A. putredinis population and finds
# ~20% meet the conventional same-strain threshold (popANI >= 0.99999). That
# number is only informative relative to a background: if two UNRELATED
# people's A. putredinis populations are typically just as similar (because
# the species has limited within-clade diversity), 20% "same strain" would
# not indicate real persistence at all. This job builds that background.
#
# Design: one reference ("anchor") MAG per named clade — the best-quality bin
# in that clade, picked by 3_make_between_person_manifest.R — against which
# every other participant IN THE SAME CLADE is compared, using their own
# baseline reads only (a "star" design: anchor vs N others, not all pairs,
# which would be ~3160 comparisons for Clade I alone). Both the anchor's own
# baseline reads and each other participant's baseline reads are mapped to
# the anchor genome and profiled independently, then compared with inStrain,
# exactly as in 2_run_instrain_compare.sh but with a FIXED reference genome
# per clade instead of each participant's own MAG.
#
# Once this job's output is collected, compare its popANI/genetic-distance
# distribution (within-clade, between-person) against the within-person
# distribution from 1_instrain_strain_retention.R — e.g. Wilcoxon rank-sum on
# (1 - popANI). If within-person distances are significantly smaller, that
# supports real persistence/microevolution for the ~20% "same strain" cases,
# and situates the rest on a genuine same-clade divergence scale rather than
# a full-strain-replacement one.
#
# NOTE on efficiency: because each task is an isolated SLURM array element,
# the anchor's own baseline reads are re-mapped and re-profiled independently
# in every task for that clade (~40 times for Clade I) rather than once and
# reused. This wastes CPU-hours but keeps the job simple and avoids any
# race condition around a shared, concurrently-written anchor profile. If
# turnaround time becomes a problem, the fix is a two-stage pipeline (profile
# every clade's anchor once in a small 5-task array, then a second array that
# only profiles the "other" participant and reuses the saved anchor profile).
#
# One array task = one (clade, other participant) row. Reads the manifest
# produced by 3_make_between_person_manifest.R (130 rows: up to 40 per clade).
#
# Container: same requirements as 2_run_instrain_compare.sh — set CONTAINER
# below (or export it) to a .sif providing inStrain, bowtie2 and samtools.
#
# Submit from the project root on Snellius:
#   sbatch scripts/4_run_instrain_between_person.sh
# Check the paths and the container first with:
#   PREFLIGHT=1 bash scripts/4_run_instrain_between_person.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG — verify these against the layout on Snellius before submitting.
# Mirrors 2_run_instrain_compare.sh; update both together if paths change.
# ---------------------------------------------------------------------------
BASE_DIR="/projects/0/prjs1784/heliuspaired"
MANIFEST="${BASE_DIR}/helius_longitudinal/results/3_species_change/5_strain_stability/instrain_between_person_manifest.csv"
OUT_DIR="${BASE_DIR}/instrain_ap_between"   # separate from instrain_ap/ (within-person)
                                             # to avoid subject_id filename collisions

# Per-batch nf-core/mag result directories holding the refined bin FASTAs.
# (A function rather than an associative array, so this also runs under the
# bash 3.2 on macOS when preflighting locally.)
bin_dir_for_batch() {
  case "$1" in
    1) echo "${BASE_DIR}/mag_results/results_batch1/GenomeBinning/DASTool/bins" ;;
    2) echo "${BASE_DIR}/mag_results/results_batch2/GenomeBinning/DASTool/bins" ;;
    3) echo "${BASE_DIR}/mag_results/results_batch3/GenomeBinning/DASTool/bins" ;;
    *) die "unknown batch: $1" ;;
  esac
}

# Where the QC'd (host-removed) reads live, and how they are named.
# Candidate patterns are tried in order; {SAMPLE}, {ID}, {TP} and {R} are
# substituted ({SAMPLE}=HELIBA_100043 -> {ID}=100043, {TP}=BA; {R} is 1 or 2).
# The first pattern that matches both mates is used. nf-core/mag did not keep
# the host-removed fastqs, so in practice the raw reads are used; human reads
# do not align to a bacterial MAG at MIN_READ_ANI, so this is harmless.
READS_DIRS=(
  "${BASE_DIR}/mag_results/results_batch1/QC_shortreads/remove_host"
  "${BASE_DIR}/mag_results/results_batch2/QC_shortreads/remove_host"
  "${BASE_DIR}/mag_results/results_batch3/QC_shortreads/remove_host"
  "${BASE_DIR}/data/fastqsraw"
)
READS_PATTERNS=(
  "HELIUS.Metagenome.{ID}.Fecal.{TP}.NA.{R}.fq.gz"
  "{SAMPLE}_run0_host_removed_{R}.fastq.gz"
  "{SAMPLE}_run1_host_removed_{R}.fastq.gz"
  "{SAMPLE}_host_removed_{R}.fastq.gz"
  "{SAMPLE}_run0_phix_removed_{R}.fastq.gz"
)

# Same thresholds as 2_run_instrain_compare.sh, so the two popANI
# distributions are directly comparable.
MIN_READ_ANI=0.95
MIN_COV=5
MIN_BREADTH=0.5

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Container image. It must provide inStrain, bowtie2 and samtools; the
# preflight below checks all three and names any that are missing.
# Override without editing this file by exporting CONTAINER=/path/to/image.sif
CONTAINER="${CONTAINER:-/home/bverhaar/singularity_images/instrain_bowtie2_samtools.img}"
# On Snellius you may need to load the module in your submit environment:
#   module load 2023 && module load Apptainer/1.2.5-GCCcore-12.3.0
APPTAINER_BIN="${APPTAINER_BIN:-$(command -v apptainer || command -v singularity || true)}"

# ---------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 1; }

# Run a tool inside the container. BIND_DIRS is extended once WORK exists.
# --no-mount hostfs: Snellius mounts every host filesystem into containers,
# including /opt, which hides the image's /opt/conda (where the tools live).
BIND_DIRS="${BASE_DIR}"
sif_exec() {
  "$APPTAINER_BIN" exec --cleanenv --no-mount hostfs -B "$BIND_DIRS" "$CONTAINER" "$@"
}

# Find the read pair for a sample; echoes "R1 R2" or dies with the paths tried.
find_reads() {
  local sample="$1" dir pat base r1 r2 tried=""
  local id="${sample#*_}" tp="${sample:4:2}"
  for dir in "${READS_DIRS[@]}"; do
    for pat in "${READS_PATTERNS[@]}"; do
      base="${pat//\{SAMPLE\}/$sample}"; base="${base//\{ID\}/$id}"; base="${base//\{TP\}/$tp}"
      r1="${dir}/${base//\{R\}/1}"
      r2="${dir}/${base//\{R\}/2}"
      if [[ -s "$r1" && -s "$r2" ]]; then echo "$r1 $r2"; return 0; fi
      tried+="  ${r1}"$'\n'
    done
  done
  die "no read pair found for ${sample}. Tried:"$'\n'"${tried}"
}

[[ -s "$MANIFEST" ]] || die "manifest not found: ${MANIFEST} (run 3_make_between_person_manifest.R and copy it here)"

# ---- Preflight: check row 1 resolves, print what was found, then stop ----
if [[ "${PREFLIGHT:-0}" == "1" ]]; then
  echo "Manifest: ${MANIFEST} ($(( $(wc -l < "$MANIFEST") - 1 )) comparisons)"
  row=$(sed -n '2p' "$MANIFEST" | tr -d '"')
  IFS=, read -r clade a_subj a_bin a_batch a_sample a_compl a_cont o_subj o_sample <<< "$row"
  echo "First row: ${clade} — anchor ${a_subj} (bin ${a_bin}, batch ${a_batch}) vs other ${o_subj}"
  fasta="$(bin_dir_for_batch "$a_batch")/${a_bin}.fa"
  [[ -s "$fasta" ]] && echo "  anchor MAG FASTA   OK  ${fasta}" \
                    || echo "  anchor MAG FASTA   MISSING  ${fasta}"
  for s in "$a_sample" "$o_sample"; do
    if reads=$(find_reads "$s" 2>/dev/null); then
      echo "  reads ${s}  OK  ${reads%% *}"
    else
      echo "  reads ${s}  MISSING — adjust READS_DIRS / READS_PATTERNS"
    fi
  done
  if [[ -z "$APPTAINER_BIN" ]]; then
    echo "  apptainer/singularity NOT on PATH — load the module first, e.g."
    echo "     module load 2023 && module load Apptainer/1.2.5-GCCcore-12.3.0"
  else
    echo "  container runtime  OK  ${APPTAINER_BIN}"
    if [[ -s "$CONTAINER" ]]; then
      echo "  image  OK  ${CONTAINER}"
      for tool in inStrain bowtie2 bowtie2-build samtools; do
        if sif_exec sh -c "command -v $tool" >/dev/null 2>&1; then
          echo "    ${tool}  OK  $(sif_exec "$tool" --version 2>&1 | head -1)"
        else
          echo "    ${tool}  MISSING from the image"
        fi
      done
    else
      echo "  image  MISSING  ${CONTAINER} — set CONTAINER=/path/to/image.sif"
    fi
  fi
  exit 0
fi

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "not a job array; submit with sbatch"

# Fail here rather than halfway through mapping if the runtime or image is missing
[[ -n "$APPTAINER_BIN" ]] || die "apptainer/singularity not on PATH — load the module before submitting"
[[ -s "$CONTAINER" ]] || die "container image not found: ${CONTAINER} (set CONTAINER=/path/to/image.sif)"
for tool in inStrain bowtie2 bowtie2-build samtools; do
  sif_exec sh -c "command -v $tool" >/dev/null 2>&1 \
    || die "'${tool}' not found inside ${CONTAINER}"
done

# ---- Manifest row for this array task (row 1 = header) ----
row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | tr -d '"')
[[ -n "$row" ]] || die "no manifest row for array index ${SLURM_ARRAY_TASK_ID}"
IFS=, read -r CLADE ANCHOR_SUBJECT ANCHOR_BIN ANCHOR_BATCH ANCHOR_SAMPLE \
              ANCHOR_COMPL ANCHOR_CONT OTHER_SUBJECT OTHER_SAMPLE <<< "$row"

echo "=== ${CLADE}: anchor ${ANCHOR_SUBJECT} (bin ${ANCHOR_BIN}, batch ${ANCHOR_BATCH}) vs other ${OTHER_SUBJECT}"
echo "    anchor sample ${ANCHOR_SAMPLE}  other sample ${OTHER_SAMPLE}"

# Named by OTHER_SUBJECT: unique per row, since each participant belongs to
# exactly one clade and appears as "other" at most once across the manifest.
RESULT="${OUT_DIR}/compare/${OTHER_SUBJECT}_genomeWide_compare.tsv"
if [[ -s "$RESULT" ]]; then
  echo "already done, skipping: ${RESULT}"
  exit 0
fi

FASTA="$(bin_dir_for_batch "$ANCHOR_BATCH")/${ANCHOR_BIN}.fa"
[[ -s "$FASTA" ]] || die "anchor MAG FASTA not found: ${FASTA}"

read -r R1_A R2_A <<< "$(find_reads "$ANCHOR_SAMPLE")"
read -r R1_O R2_O <<< "$(find_reads "$OTHER_SAMPLE")"

# Work on node-local scratch; only the results are copied back.
WORK="${TMPDIR:-/scratch-local}/instrain_btw_${OTHER_SUBJECT}_$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "${OUT_DIR}/compare" "${OUT_DIR}/profile"

# The container needs the project tree (bins, reads) and the node-local work
# dir. The read files may sit outside BASE_DIR, so bind their parents too.
BIND_DIRS="${BASE_DIR},${WORK},$(dirname "$R1_A"),$(dirname "$R1_O")"

# ---- Scaffold-to-bin file: every scaffold belongs to the anchor genome ----
STB="${WORK}/${ANCHOR_BIN}.stb"
awk -v g="${ANCHOR_BIN}.fa" '/^>/{print substr($1,2)"\t"g}' "$FASTA" > "$STB"
echo "    scaffolds in anchor MAG: $(wc -l < "$STB")"

# ---- Map anchor's own reads AND the other participant's reads to the SAME
#      anchor genome (unlike 2_run_instrain_compare.sh, the reference here is
#      fixed per clade, not the read-owner's own MAG) ----
IDX="${WORK}/idx"
sif_exec bowtie2-build --threads "$THREADS" -q "$FASTA" "$IDX"

map_sample() {
  local tag="$1" r1="$2" r2="$3" bam="${WORK}/${1}.bam"
  sif_exec bowtie2 -x "$IDX" -1 "$r1" -2 "$r2" -p "$THREADS" \
      2> "${WORK}/${tag}.bowtie2.log" \
    | sif_exec samtools sort -@ "$THREADS" -o "$bam" -
  sif_exec samtools index "$bam"
  echo "    ${tag} alignment rate: $(grep 'overall alignment rate' "${WORK}/${tag}.bowtie2.log" || echo NA)"
}
map_sample anchor "$R1_A" "$R2_A"
map_sample other  "$R1_O" "$R2_O"

# ---- inStrain profile, one per sample against the anchor genome ----
for tag in anchor other; do
  sif_exec inStrain profile "${WORK}/${tag}.bam" "$FASTA" \
    -o "${WORK}/${tag}.IS" \
    -p "$THREADS" \
    -s "$STB" \
    --min_read_ani "$MIN_READ_ANI" \
    --skip_plot_generation
done

# ---- inStrain compare: anchor vs other participant ----
sif_exec inStrain compare \
  -i "${WORK}/anchor.IS" "${WORK}/other.IS" \
  -o "${WORK}/compare.IS" \
  -p "$THREADS" \
  -s "$STB" \
  --min_cov "$MIN_COV" \
  --breadth "$MIN_BREADTH" \
  --skip_plot_generation

# ---- Collect: prepend clade + both subject ids so the pooled table is
#      self-documenting without re-parsing filenames ----
GW=$(find "${WORK}/compare.IS/output" -name "*genomeWide_compare.tsv" | head -1)
[[ -s "$GW" ]] || die "inStrain compare produced no genomeWide table for ${OTHER_SUBJECT}"
awk -v c="$CLADE" -v a="$ANCHOR_SUBJECT" -v o="$OTHER_SUBJECT" \
  'NR==1{print "clade\tanchor_subject_id\tother_subject_id\t"$0; next}
   {print c"\t"a"\t"o"\t"$0}' "$GW" > "$RESULT"

CMP=$(find "${WORK}/compare.IS/output" -name "*comparisonsTable.tsv" | head -1)
[[ -s "$CMP" ]] && awk -v c="$CLADE" -v a="$ANCHOR_SUBJECT" -v o="$OTHER_SUBJECT" \
  'NR==1{print "clade\tanchor_subject_id\tother_subject_id\t"$0; next}
   {print c"\t"a"\t"o"\t"$0}' "$CMP" \
  > "${OUT_DIR}/compare/${OTHER_SUBJECT}_comparisonsTable.tsv"

# Keep the per-sample genome_info for coverage/breadth diagnostics
for tag in anchor other; do
  GI=$(find "${WORK}/${tag}.IS/output" -name "*genome_info.tsv" | head -1)
  [[ -s "$GI" ]] && awk -v o="$OTHER_SUBJECT" -v t="$tag" \
    'NR==1{print "other_subject_id\trole\t"$0; next}{print o"\t"t"\t"$0}' "$GI" \
    > "${OUT_DIR}/profile/${OTHER_SUBJECT}_${tag}_genome_info.tsv"
done

echo "=== done: ${RESULT}"
