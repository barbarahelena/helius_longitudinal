#!/usr/bin/env python3
"""
Convert bin_summary.tsv from wide format to long format.
Extracts bin IDs and depth columns, reshaping to long format.
Can combine multiple bin_summary files from different batches.
"""

import pandas as pd
import sys

def convert_to_long_format(input_files, output_file):
    """
    Convert wide format bin summary to long format with columns:
    bin_id, sample, depth
    
    Args:
        input_files: List of paths to input bin_summary.tsv files
        output_file: Path to output CSV file
    """
    all_dfs = []
    
    for input_file in input_files:
        print(f"\nProcessing {input_file}...")
        
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
    
        # Get the bin column (first column with bin names)
        bin_col = 'bin'
        
        # Find all depth columns (columns starting with 'Depth ')
        depth_cols = [col for col in df.columns if col.startswith('Depth ')]
        
        # Extract sample names from depth column names
        # 'Depth HELIBA_100245' -> 'HELIBA_100245'
        sample_names = [col.replace('Depth ', '') for col in depth_cols]
        
        # Select only bin and depth columns
        df_subset = df[[bin_col] + depth_cols].copy()
        
        # Filter out unbinned bins (they are low quality leftovers)
        n_before = len(df_subset)
        df_subset = df_subset[~df_subset[bin_col].str.contains('Unbinned', case=False, na=False)]
        n_filtered = n_before - len(df_subset)
        if n_filtered > 0:
            print(f"  Filtered out {n_filtered} unbinned bin(s)")
        
        # Melt/reshape to long format
        df_long = df_subset.melt(
            id_vars=[bin_col],
            value_vars=depth_cols,
            var_name='sample',
            value_name='depth'
        )
        
        # Clean up sample column (remove 'Depth ' prefix)
        df_long['sample'] = df_long['sample'].str.replace('Depth ', '')
        
        # Rename bin column to bin_id
        df_long = df_long.rename(columns={bin_col: 'bin_id'})
        
        # Strip .fa extension from bin_id
        df_long['bin_id'] = df_long['bin_id'].str.replace('.fa$', '', regex=True)
        
        # Remove rows with missing depth values
        df_long = df_long.dropna(subset=['depth'])
        
        print(f"  Found {len(df_subset)} bins with {len(depth_cols)} samples")
        print(f"  Generated {len(df_long)} depth entries")
        
        all_dfs.append(df_long)
    
    # Combine all dataframes
    print(f"\nCombining data from {len(input_files)} file(s)...")
    df_combined = pd.concat(all_dfs, ignore_index=True)
    
    # Remove any duplicate entries (same bin_id, sample combination)
    n_before_dedup = len(df_combined)
    df_combined = df_combined.drop_duplicates(subset=['bin_id', 'sample'], keep='first')
    n_duplicates = n_before_dedup - len(df_combined)
    if n_duplicates > 0:
        print(f"  Removed {n_duplicates} duplicate entries")
    
    # Reorder columns
    df_combined = df_combined[['bin_id', 'sample', 'depth']]
    
    # Sort by bin_id and sample
    df_combined = df_combined.sort_values(['bin_id', 'sample']).reset_index(drop=True)
    
    # Save to CSV
    df_combined.to_csv(output_file, index=False)
    
    print(f"\nOutput: {len(df_combined)} rows written to {output_file}")
    print(f"\nFirst few rows:")
    print(df_combined.head(10))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert_bin_summary_to_long.py <input_file1> [input_file2 ...] <output_file>")
        print("\nArguments:")
        print("  input_file(s) : One or more bin_summary.tsv files to process")
        print("  output_file   : Path to output CSV file (last argument)")
        print("\nExamples:")
        print("  Single file:")
        print("    python convert_bin_summary_to_long.py bin_summary.tsv bin_depths_long.csv")
        print("\n  Multiple batches:")
        print("    python convert_bin_summary_to_long.py batch1/bin_summary.tsv batch2/bin_summary.tsv batch3/bin_summary.tsv bin_depths_long.csv")
        sys.exit(1)
    
    # All arguments except the last are input files
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    
    print(f"Processing {len(input_files)} input file(s):")
    for f in input_files:
        print(f"  - {f}")
    
    convert_to_long_format(input_files, output_file)
