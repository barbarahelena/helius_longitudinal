#!/usr/bin/env python3
"""
Create input sheet for bins from bin_summary.tsv.
Generates a CSV with: bin_id, bin_fasta, checkm_completeness, checkm_contamination
Can combine multiple bin_summary files from different batches.
"""

import pandas as pd
import sys
import os

def create_bin_input_sheet(input_files, output_file, bins_dirs=None):
    """
    Create input sheet from bin_summary.tsv
    
    Args:
        input_files: List of paths to input bin_summary.tsv files
        output_file: Path to output CSV file
        bins_dirs: List of directories containing the bin FASTA files (optional, will be inferred if not provided)
    """
    all_dfs = []
    
    for i, input_file in enumerate(input_files):
        print(f"\nProcessing {input_file}...")
        
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
    
        # Select relevant columns
        # Column names: bin, Name, Completeness, Contamination
        required_cols = ['bin', 'Name', 'Completeness', 'Contamination']
        
        # Check if all required columns exist
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: Missing columns: {missing_cols}")
            print(f"Available columns: {df.columns.tolist()}")
            sys.exit(1)
        
        # Create output dataframe
        df_output = df[required_cols].copy()
        
        # Rename columns to match output format
        df_output = df_output.rename(columns={
            'bin': 'bin_id',
            'Name': 'bin_name',
            'Completeness': 'checkm_completeness',
            'Contamination': 'checkm_contamination'
        })
        
        # Strip .fa extension from bin_id
        df_output['bin_id'] = df_output['bin_id'].str.replace('.fa$', '', regex=True)
        
        # Filter out unbinned bins (they are low quality leftovers)
        n_before = len(df_output)
        df_output = df_output[~df_output['bin_name'].str.contains('Unbinned', case=False, na=False)]
        n_filtered = n_before - len(df_output)
        if n_filtered > 0:
            print(f"  Filtered out {n_filtered} unbinned bin(s)")
        
        # Infer bins directory if not provided
        if bins_dirs is None:
            # Get directory of input file
            input_dir = os.path.dirname(os.path.abspath(input_file))
            # Assume bins are in DASTool/bins subdirectory
            bins_dir = os.path.join(input_dir, 'DASTool', 'bins')
            print(f"  Using bins directory: {bins_dir}")
        else:
            bins_dir = bins_dirs[i] if i < len(bins_dirs) else bins_dirs[0]
            print(f"  Using bins directory: {bins_dir}")
        
        # Create bin_fasta paths
        # The bin name (e.g., MEGAHIT-CONCOCTRefined-group-10.13) should match the fasta filename
        df_output['bin_fasta'] = df_output['bin_name'].apply(
            lambda x: os.path.join(bins_dir, f"{x}.fa")
        )
        
        # Reorder columns
        df_output = df_output[['bin_id', 'bin_fasta', 'checkm_completeness', 'checkm_contamination']]
        
        print(f"  Found {len(df_output)} bins")
        
        all_dfs.append(df_output)
    
    # Combine all dataframes
    print(f"\nCombining data from {len(input_files)} file(s)...")
    df_combined = pd.concat(all_dfs, ignore_index=True)
    
    # Remove any duplicate bin_ids (keep first occurrence)
    n_before_dedup = len(df_combined)
    df_combined = df_combined.drop_duplicates(subset=['bin_id'], keep='first')
    n_duplicates = n_before_dedup - len(df_combined)
    if n_duplicates > 0:
        print(f"  Warning: Removed {n_duplicates} duplicate bin_id(s)")
    
    # Sort by bin_id
    df_combined = df_combined.sort_values('bin_id').reset_index(drop=True)
    
    # Save to CSV
    df_combined.to_csv(output_file, index=False)
    
    print(f"\nCreated input sheet with {len(df_combined)} bins")
    print(f"Output written to: {output_file}")
    print(f"\nFirst few rows:")
    print(df_combined.head(10))
    
    # Summary statistics
    print(f"\nSummary:")
    print(f"  Completeness: mean={df_combined['checkm_completeness'].mean():.2f}, "
          f"median={df_combined['checkm_completeness'].median():.2f}")
    print(f"  Contamination: mean={df_combined['checkm_contamination'].mean():.2f}, "
          f"median={df_combined['checkm_contamination'].median():.2f}")
    
    # Count high-quality bins (>90% complete, <5% contamination)
    high_quality = df_combined[
        (df_combined['checkm_completeness'] >= 90) & 
        (df_combined['checkm_contamination'] < 5)
    ]
    print(f"  High-quality bins (≥90% complete, <5% contamination): {len(high_quality)}")
    
    # Medium-quality bins (>50% complete, <10% contamination)
    medium_quality = df_combined[
        (df_combined['checkm_completeness'] >= 50) & 
        (df_combined['checkm_contamination'] < 10)
    ]
    print(f"  Medium-quality bins (≥50% complete, <10% contamination): {len(medium_quality)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python create_bin_input_sheet.py <input_file1> [input_file2 ...] <output_file> [bins_dir1 bins_dir2 ...]")
        print("\nArguments:")
        print("  input_file(s) : One or more bin_summary.tsv files to process")
        print("  output_file   : Path to output CSV file")
        print("  bins_dir(s)   : (Optional) Directory(ies) containing bin FASTA files")
        print("                  If not provided, assumes DASTool/bins relative to each input file")
        print("\nExamples:")
        print("  Single file:")
        print("    python create_bin_input_sheet.py bin_summary.tsv bin_input_sheet.csv")
        print("\n  Multiple batches (auto-detect bins directories):")
        print("    python create_bin_input_sheet.py batch1/bin_summary.tsv batch2/bin_summary.tsv batch3/bin_summary.tsv bin_input_sheet.csv")
        print("\n  Multiple batches (specify bins directories):")
        print("    python create_bin_input_sheet.py batch1/bin_summary.tsv batch2/bin_summary.tsv bin_input_sheet.csv batch1/bins batch2/bins")
        sys.exit(1)
    
    # Determine if bins directories are provided
    # Last argument is always output file
    # If there are more args than input files + 1, the extra ones are bins_dirs
    output_file = None
    input_files = []
    bins_dirs = None
    
    # Simple heuristic: if an argument ends with .csv, it's likely the output file
    csv_indices = [i for i, arg in enumerate(sys.argv[1:], 1) if arg.endswith('.csv')]
    
    if len(csv_indices) == 1:
        output_idx = csv_indices[0]
        input_files = sys.argv[1:output_idx]
        output_file = sys.argv[output_idx]
        if len(sys.argv) > output_idx + 1:
            bins_dirs = sys.argv[output_idx + 1:]
    else:
        # Fallback: assume last argument is output, rest are inputs
        input_files = sys.argv[1:-1]
        output_file = sys.argv[-1]
    
    print(f"Processing {len(input_files)} input file(s):")
    for f in input_files:
        print(f"  - {f}")
    
    if bins_dirs:
        print(f"\nUsing {len(bins_dirs)} bins directory(ies):")
        for d in bins_dirs:
            print(f"  - {d}")
    
    create_bin_input_sheet(input_files, output_file, bins_dirs)
