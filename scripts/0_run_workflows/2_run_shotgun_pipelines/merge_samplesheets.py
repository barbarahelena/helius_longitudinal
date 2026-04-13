#!/usr/bin/env python3
"""
Merge multiple samplesheets and remove the last two columns.
Combines samplesheet.csv, samplesheet_2.csv, and samplesheet_3.csv
"""

import pandas as pd
import sys
import os

def merge_samplesheets(input_files, output_file):
    """
    Merge multiple samplesheet CSV files and remove last two columns
    
    Args:
        input_files: List of paths to input samplesheet CSV files
        output_file: Path to output merged CSV file
    """
    all_dfs = []
    
    for input_file in input_files:
        print(f"Reading {input_file}...")
        
        # Check if file exists
        if not os.path.exists(input_file):
            print(f"  Warning: File not found, skipping: {input_file}")
            continue
        
        # Read the CSV file
        df = pd.read_csv(input_file)
        print(f"  Found {len(df)} rows with {len(df.columns)} columns")
        print(f"  Columns: {', '.join(df.columns.tolist())}")
        
        all_dfs.append(df)
    
    if not all_dfs:
        print("Error: No valid input files found!")
        sys.exit(1)
    
    # Combine all dataframes
    print(f"\nCombining data from {len(all_dfs)} file(s)...")
    df_combined = pd.concat(all_dfs, ignore_index=True)
    print(f"  Total rows before deduplication: {len(df_combined)}")
    
    # Remove duplicate rows (based on all columns)
    n_before_dedup = len(df_combined)
    df_combined = df_combined.drop_duplicates()
    n_duplicates = n_before_dedup - len(df_combined)
    if n_duplicates > 0:
        print(f"  Removed {n_duplicates} duplicate row(s)")
    
    # Remove last two columns
    original_cols = df_combined.columns.tolist()
    if len(original_cols) >= 2:
        df_combined = df_combined.iloc[:, :-2]
        removed_cols = original_cols[-2:]
        print(f"\nRemoved last two columns: {', '.join(removed_cols)}")
        print(f"Remaining columns: {', '.join(df_combined.columns.tolist())}")
    else:
        print(f"\nWarning: DataFrame has only {len(original_cols)} column(s), cannot remove last two columns")
    
    # Sort by sample name
    if 'sample' in df_combined.columns:
        df_combined = df_combined.sort_values('sample').reset_index(drop=True)
        print(f"Sorted by 'sample' column")
    
    # Save to CSV
    df_combined.to_csv(output_file, index=False)
    
    print(f"\nMerged samplesheet with {len(df_combined)} rows written to: {output_file}")
    print(f"\nFirst few rows:")
    print(df_combined.head(10))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_samplesheets.py <input_file1> [input_file2 ...] <output_file>")
        print("\nArguments:")
        print("  input_file(s) : One or more samplesheet CSV files to merge")
        print("  output_file   : Path to output merged CSV file (last argument)")
        print("\nExample:")
        print("  python merge_samplesheets.py samplesheet.csv samplesheet_2.csv samplesheet_3.csv merged_samplesheet.csv")
        sys.exit(1)
    
    # All arguments except the last are input files
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    
    print(f"Input files to merge:")
    for f in input_files:
        print(f"  - {f}")
    print(f"Output file: {output_file}\n")
    
    merge_samplesheets(input_files, output_file)
