#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --job-name=alistipes_phylo_c70
#SBATCH --output=alistipes_phylo_c70_%j.out

# Alistipes phylogenomics (Panaroo -> IQ-TREE3) on bins with CheckM2
# completeness >= 70% (contamination < 10%), reusing the existing Bakta GFF3s.
# Prepare the samplesheet first:  python3 scripts/make_alistipes_samplesheet_c70.py
# Results: phylogenomics/alistipes_c70/{panaroo,iqtree,pipeline_info}

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run /projects/0/prjs1784/heliuspaired/workflows/phylo/main.nf \
    -profile snellius \
    --taxon       alistipes_c70 \
    --samplesheet data/alistipes_bakta_samplesheet_c70.csv \
    --bakta_input true \
    --core_thresh 0.98 \
    --bootstrap   1000 \
    -resume
