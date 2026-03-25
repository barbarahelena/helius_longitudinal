#!/bin/bash
#SBATCH --cpus-per-task=2       # Number of cores
#SBATCH --mem=16G                # Memory per node
#SBATCH --time=20:00:00         # Runtime limit (HH:MM:SS)
#SBATCH --job-name=metabo       # Job name
#SBATCH --output=metabo_%j.out  # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run pipelines/metaboflow/main.nf \
    --input bin_input_sheet.csv \
    --depths bin_depths_long.csv \
    -profile snellius \
    --dram_db databases/dram_db \
    --outdir metabo_results \
    --email bjh.verhaar@gmail.com \
    -resume