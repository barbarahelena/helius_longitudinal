#!/bin/bash
python /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/scripts/3_species_change/5_mlmodels/2_ml_model/XGBeast_adj.py \
    -name eth_baseline \
    -path /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/results/3_species_change/5_mlmodels/eth_base \
    -x class \
    -n 5 \
    -rand_seed 4321 \
    -param /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/scripts/3_species_change/5_mlmodels/2_ml_model/param_grid_mb.json
python /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/scripts/3_species_change/5_mlmodels/2_ml_model/XGBeast_adj.py \
    -name eth_followup \
    -path /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/results/3_species_change/5_mlmodels/eth_fu \
    -x class \
    -n 5 \
    -rand_seed 4321 \
    -param /Users/barbaraverhaar/Documents/helius_longitudinal/helius_longitudinal/scripts/3_species_change/5_mlmodels/2_ml_model/param_grid_mb.json
