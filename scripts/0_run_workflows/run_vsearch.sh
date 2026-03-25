#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate nextflow
nextflow run pipelines/vsearchpipeline/main.nf \
        -profile snellius \
        --input data/samplesheet_combined.csv \
        --outdir helius16s/results \
        --cluster_minsize 10 \
        --cluster_alpha 2.0 \
        --run_decontam true \
        --rarelevel 13000 \
        -resume