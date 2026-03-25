#!/bin/bash
#SBATCH --cpus-per-task=4      # Number of cores
#SBATCH --mem=8G                 # Memory per node
#SBATCH --time=72:00:00          # Runtime limit (HH:MM:SS)
#SBATCH --job-name=mag        # Job name
#SBATCH --output=mag_%j.out     # Standard output file

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run pipelines/mag/main.nf \
    --input data/samplesheet_3.csv \
    -profile snellius \
    --host_genome GRCh38 \
    --coassemble_group \
    --bbnorm true \
    --bbnorm_target 40 \
    --bbnorm_min 2 \
    --skip_comebin \
    --skip_spades \
    --skip_metaeuk \
    --skip_prokka \
    --skip_prodigal \
    --binqc_tool checkm2 \
    --refine_bins_dastool \
    --postbinning_input refined_bins_only \
    --gtdb_db databases/release226 \
    --checkm2_db databases/CheckM2_database/uniref100.KO.1.dmnd \
    --outdir results_batch3 \
    --email bjh.verhaar@gmail.com \
    -resume