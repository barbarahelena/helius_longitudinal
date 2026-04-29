#!/bin/bash
#SBATCH --cpus-per-task=2       # Number of cores
#SBATCH --mem=8G                # Memory per node
#SBATCH --time=6:00:00         # Runtime limit (HH:MM:SS)
#SBATCH --job-name=senegal       # Job name
#SBATCH --output=senegal_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run /projects/prjs1784/pipelines/annotationpipeline/main.nf \
    --input data/senegalimassilia_samplesheet.csv \
    -profile snellius \
    --bakta_database /projects/prjs1784/heliuspaired/annotation_parabacteroides/bakta/db \
    --outdir annotation_senegalimassilia \
    --annotation true \
    --email bjh.verhaar@gmail.com \
    -resume