import pandas as pd
import h5py
import numpy as np
import sys

# --- Check if the correct number of arguments is provided ---
if len(sys.argv) != 4:
    print("Usage: python script_name.py <path_to_rna_tsv> <path_to_atac_tsv> <output_h5_file>")
    sys.exit(1)

# --- Get the file paths and output filename from the command line arguments ---
rna_tsv_file = sys.argv[1]
atac_tsv_file = sys.argv[2]
output_h5_file = sys.argv[3]

# --- Function to read a TSV file into a Pandas DataFrame ---
def read_tsv(file_path):
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    return df

# --- Read the TSV files ---
try:
    rna_df = read_tsv(rna_tsv_file)
    atac_df = read_tsv(atac_tsv_file)
except FileNotFoundError:
    print("Error: One or both of the specified files were not found.")
    sys.exit(1)

# --- Create the H5 file ---
with h5py.File(output_h5_file, 'w') as hf:
    # --- Store scRNA-seq data ---
    rna_group = hf.create_group('rna')
    rna_group.create_dataset('barcodes', data=np.array(rna_df.columns, dtype='S'))
    rna_group.create_dataset('features', data=np.array(rna_df.index, dtype='S'))
    rna_group.create_dataset('data', data=rna_df.values.T)

    # --- Store scATAC-seq data ---
    atac_group = hf.create_group('atac')
    atac_group.create_dataset('barcodes', data=np.array(atac_df.columns, dtype='S'))
    atac_group.create_dataset('features', data=np.array(atac_df.index, dtype='S'))
    atac_group.create_dataset('data', data=atac_df.values.T)

    # --- (Optional) Store joint information or metadata if available ---
    joint_group = hf.create_group('joint')
    # ...

print(f"Successfully created joint H5 file: {output_h5_file}")