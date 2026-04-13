#!/bin/bash
#SBATCH --cpus-per-task=2          # Number of cores
#SBATCH --mem=8G                   # Memory per node
#SBATCH --time=72:00:00            # Runtime limit (HH:MM:SS)
#SBATCH --job-name=cayman          # Job name
#SBATCH --output=cayman_%j.out     # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run /projects/prjs1784/pipelines/caymanflow \
  --input /projects/prjs1784/heliuspaired/data/cayman_samplesheet.csv\
  --outdir caymanflow/caymanflow_results \
  --genome GRCh38 \
  --bwa_index /projects/prjs1784/heliuspaired/caymanflow/bwa_index \
  --cayman_database /projects/prjs1784/heliuspaired/caymanflow/cayman_unzip/human-gut.fna.gz \
  --cayman_annotations /projects/prjs1784/heliuspaired/caymanflow/cayman_unzip/human-gut_annotations.csv \
  -work-dir /scratch-shared/bverhaar/work_caymanflow \
  -profile snellius \
  -resume