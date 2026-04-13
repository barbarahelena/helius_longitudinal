#!/bin/bash
#SBATCH --cpus-per-task=2       # Number of cores
#SBATCH --mem=8G                # Memory per node
#SBATCH --time=20:00:00         # Runtime limit (HH:MM:SS)
#SBATCH --job-name=anno       # Job name
#SBATCH --output=anno_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run pipelines/annotationpipeline/main.nf \
    --input data/alistipes_samplesheet.csv \
    -profile snellius \
    --outdir annotation_alistipes \
    --annotation true \
    --email bjh.verhaar@gmail.com \
    -resume