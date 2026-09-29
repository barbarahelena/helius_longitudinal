#!/bin/bash
#SBATCH --job-name=instrain_ap
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --array=1-126%20
#SBATCH --output=logs/instrain_%A_%a.out
#SBATCH --error=logs/instrain_%A_%a.err

# Alistipes putredinis strain retention: baseline vs follow-up, per participant.
#
# Why: MAGs were co-assembled per participant from pooled baseline + follow-up
# reads and refined with DAS Tool, which emits a non-redundant bin set. Each
# participant therefore has exactly one A. putredinis MAG, so "clade at
# baseline" and "clade at follow-up" are read off the same genome and can
# never differ. This job measures strain retention directly: it maps each
# timepoint's reads to that participant's own MAG and compares the two
# populations with inStrain (popANI).
#
# One array task = one participant. Reads the manifest produced by
# 1_make_instrain_manifest.R; rows are sorted so the adequate-depth
# participants come first, hence --array=1-126 above. Set --array=1-164 to
# include the low-coverage ones as well (they are flagged, not dropped).
#
# Container: set CONTAINER below (or export it) to a .sif providing inStrain,
# bowtie2 and samtools. Note that the plain inStrain biocontainer is
# Python-only and ships neither bowtie2 nor samtools, so an image built from
# just that will not work here — the preflight names any missing tool.
#
# Submit from the project root on Snellius:
#   sbatch scripts/2_run_instrain_compare.sh
# Check the paths and the container first with:
#   PREFLIGHT=1 bash scripts/2_run_instrain_compare.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG — verify these against the layout on Snellius before submitting
# ---------------------------------------------------------------------------
BASE_DIR="/projects/0/prjs1784/heliuspaired"
MANIFEST="${BASE_DIR}/helius_longitudinal/results/3_species_change/5_strain_stability/instrain_manifest.csv"
OUT_DIR="${BASE_DIR}/instrain_ap"

# Per-batch nf-core/mag result directories holding the refined bin FASTAs.
# The manifest's `batch` column selects which one to look in.
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

# inStrain thresholds. popANI >= 0.99999 is the conventional "same strain"
# cutoff (Olm et al. 2021); applied at the collection step, not here.
MIN_READ_ANI=0.95   # reads below this identity are not counted
MIN_COV=5           # minimum coverage for a position to enter the comparison
MIN_BREADTH=0.5     # minimum fraction of the genome covered in both samples

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

[[ -s "$MANIFEST" ]] || die "manifest not found: ${MANIFEST}"

# ---- Preflight: check row 1 resolves, print what was found, then stop ----
if [[ "${PREFLIGHT:-0}" == "1" ]]; then
  echo "Manifest: ${MANIFEST} ($(( $(wc -l < "$MANIFEST") - 1 )) participants)"
  row=$(sed -n '2p' "$MANIFEST" | tr -d '"')
  IFS=, read -r subject bin batch s_bl s_fu _ <<< "$row"
  echo "First participant: ${subject} (bin ${bin}, batch ${batch})"
  fasta="$(bin_dir_for_batch "$batch")/${bin}.fa"
  [[ -s "$fasta" ]] && echo "  MAG FASTA   OK  ${fasta}" \
                    || echo "  MAG FASTA   MISSING  ${fasta}"
  for s in "$s_bl" "$s_fu"; do
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
      # Every tool the job needs must be inside the image
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
# Fail here rather than halfway through mapping if the image is missing a tool
for tool in inStrain bowtie2 bowtie2-build samtools; do
  sif_exec sh -c "command -v $tool" >/dev/null 2>&1 \
    || die "'${tool}' not found inside ${CONTAINER}"
done

# ---- Manifest row for this array task (row 1 = header) ----
row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST" | tr -d '"')
[[ -n "$row" ]] || die "no manifest row for array index ${SLURM_ARRAY_TASK_ID}"
IFS=, read -r SUBJECT BIN BATCH SAMPLE_BL SAMPLE_FU DEPTH_BL DEPTH_FU _ <<< "$row"

echo "=== participant ${SUBJECT} | bin ${BIN} | batch ${BATCH}"
echo "    baseline ${SAMPLE_BL} (${DEPTH_BL}x)  follow-up ${SAMPLE_FU} (${DEPTH_FU}x)"

RESULT="${OUT_DIR}/compare/${SUBJECT}_genomeWide_compare.tsv"
if [[ -s "$RESULT" ]]; then
  echo "already done, skipping: ${RESULT}"
  exit 0
fi

FASTA="$(bin_dir_for_batch "$BATCH")/${BIN}.fa"
[[ -s "$FASTA" ]] || die "MAG FASTA not found: ${FASTA}"

read -r R1_BL R2_BL <<< "$(find_reads "$SAMPLE_BL")"
read -r R1_FU R2_FU <<< "$(find_reads "$SAMPLE_FU")"

# Work on node-local scratch; only the results are copied back.
WORK="${TMPDIR:-/scratch-local}/instrain_${SUBJECT}_$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "${OUT_DIR}/compare" "${OUT_DIR}/profile"

# The container needs the project tree (bins, reads) and the node-local work dir.
# The read files may sit outside BASE_DIR, so bind their parents too.
BIND_DIRS="${BASE_DIR},${WORK},$(dirname "$R1_BL"),$(dirname "$R1_FU")"

# ---- Scaffold-to-bin file: every scaffold belongs to this one genome ----
STB="${WORK}/${BIN}.stb"
awk -v g="${BIN}.fa" '/^>/{print substr($1,2)"\t"g}' "$FASTA" > "$STB"
echo "    scaffolds in MAG: $(wc -l < "$STB")"

# ---- Map each timepoint's reads to this participant's own MAG ----
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
map_sample baseline  "$R1_BL" "$R2_BL"
map_sample followup  "$R1_FU" "$R2_FU"

# ---- inStrain profile per timepoint ----
for tag in baseline followup; do
  sif_exec inStrain profile "${WORK}/${tag}.bam" "$FASTA" \
    -o "${WORK}/${tag}.IS" \
    -p "$THREADS" \
    -s "$STB" \
    --min_read_ani "$MIN_READ_ANI" \
    --skip_plot_generation
done

# ---- inStrain compare: baseline vs follow-up ----
sif_exec inStrain compare \
  -i "${WORK}/baseline.IS" "${WORK}/followup.IS" \
  -o "${WORK}/compare.IS" \
  -p "$THREADS" \
  -s "$STB" \
  --min_cov "$MIN_COV" \
  --breadth "$MIN_BREADTH" \
  --skip_plot_generation

# ---- Collect: prepend the participant id so the tables concatenate cleanly ----
GW=$(find "${WORK}/compare.IS/output" -name "*genomeWide_compare.tsv" | head -1)
[[ -s "$GW" ]] || die "inStrain compare produced no genomeWide table for ${SUBJECT}"
awk -v s="$SUBJECT" 'NR==1{print "subject_id\t"$0; next}{print s"\t"$0}' "$GW" > "$RESULT"

CMP=$(find "${WORK}/compare.IS/output" -name "*comparisonsTable.tsv" | head -1)
[[ -s "$CMP" ]] && awk -v s="$SUBJECT" 'NR==1{print "subject_id\t"$0; next}{print s"\t"$0}' "$CMP" \
  > "${OUT_DIR}/compare/${SUBJECT}_comparisonsTable.tsv"

# Keep the per-timepoint genome_info for coverage/breadth diagnostics
for tag in baseline followup; do
  GI=$(find "${WORK}/${tag}.IS/output" -name "*genome_info.tsv" | head -1)
  [[ -s "$GI" ]] && awk -v s="$SUBJECT" -v t="$tag" \
    'NR==1{print "subject_id\ttimepoint\t"$0; next}{print s"\t"t"\t"$0}' "$GI" \
    > "${OUT_DIR}/profile/${SUBJECT}_${tag}_genome_info.tsv"
done

echo "=== done: ${RESULT}"
