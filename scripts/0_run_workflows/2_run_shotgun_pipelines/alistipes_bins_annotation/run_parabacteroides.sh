#!/bin/bash
#SBATCH --cpus-per-task=2       # Number of cores
#SBATCH --mem=8G                # Memory per node
#SBATCH --time=8:00:00         # Runtime limit (HH:MM:SS)
#SBATCH --job-name=para       # Job name
#SBATCH --output=para_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run /projects/prjs1784/pipelines/annotationpipeline/main.nf \
    --input data/parabacteroides_samplesheet.csv \
    -profile snellius \
    --outdir annotation_parabacteroides \
    --annotation true \
    --email bjh.verhaar@gmail.com \
    -resume