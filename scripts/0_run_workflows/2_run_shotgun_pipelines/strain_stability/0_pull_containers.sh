#!/bin/bash
# Pull the containers for the strain-retention job once, on the login node.
#
# Run this before submitting 2_run_instrain_compare.sh. The array runs 20 tasks
# at a time and each would otherwise pull the same images from quay.io
# simultaneously, which is slow and hammers the registry.
#
#   bash 0_pull_containers.sh
#
# The inStrain biocontainer is Python-only (biopython, pysam, pandas, ...) and
# ships neither bowtie2 nor samtools, so the mapping step uses its own images
# and the pipeline is piped between them on the host.

set -euo pipefail

BASE_DIR="${BASE_DIR:-/projects/prjs1784/heliuspaired}"
SIF_DIR="${SIF_DIR:-${BASE_DIR}/containers}"

# Tags verified against quay.io/api/v1/repository/biocontainers/<tool>/tag/
INSTRAIN_URI="docker://quay.io/biocontainers/instrain:1.10.0--pyhdfd78af_0"
BOWTIE2_URI="docker://quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"
SAMTOOLS_URI="docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

# Apptainer replaced Singularity; use whichever this cluster provides.
APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
if [[ -z "$APPTAINER_BIN" ]]; then
  echo "ERROR: neither apptainer nor singularity found on PATH." >&2
  echo "       On Snellius load the module first, e.g.:" >&2
  echo "         module load 2023 && module load Apptainer/1.2.5-GCCcore-12.3.0" >&2
  echo "       (run 'module spider apptainer' to see what is available)" >&2
  exit 1
fi
echo "Using: ${APPTAINER_BIN} ($("$APPTAINER_BIN" --version))"

mkdir -p "$SIF_DIR"
# Keep the layer cache off $HOME, which is usually quota-limited
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SIF_DIR}/.cache}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_CACHEDIR"

pull() {
  local name="$1" uri="$2" sif="${SIF_DIR}/$1.sif"
  if [[ -s "$sif" ]]; then
    echo "already present: ${sif}"
  else
    echo "pulling ${name} <- ${uri}"
    "$APPTAINER_BIN" pull "$sif" "$uri"
  fi
}

pull instrain  "$INSTRAIN_URI"
pull bowtie2   "$BOWTIE2_URI"
pull samtools  "$SAMTOOLS_URI"

echo
echo "Verifying the tools run inside the images:"
"$APPTAINER_BIN" exec "${SIF_DIR}/instrain.sif" inStrain --version
"$APPTAINER_BIN" exec "${SIF_DIR}/bowtie2.sif"  bowtie2 --version | head -1
"$APPTAINER_BIN" exec "${SIF_DIR}/samtools.sif" samtools --version | head -1

echo
echo "Containers ready in: ${SIF_DIR}"
echo "Record these versions in the methods; the .sif files are the actual pin."
