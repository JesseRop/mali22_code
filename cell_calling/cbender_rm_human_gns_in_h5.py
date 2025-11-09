# This script reads the sample names from a CSV file and then iterates over each sample name to read the corresponding h5 file and save its metadata as a CSV file.
# Run as below
        ## source /software/team222/jr35/cellbender_gpu/cellbender_gpu_env/bin/activate
        ## python cbender_rm_human_gns_in_h5.py

# Load required modules and functions
from pathlib import Path
import h5py
import numpy as np
from scipy.sparse import csc_matrix
import shutil
# import os

## set working directory
# os.chdir('/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf')

def filter_h5_ensg(infile, outfile):
    with h5py.File(infile, "r") as f_in:
        # Load sparse matrix
        data = f_in["matrix"]["data"][:]
        indices = f_in["matrix"]["indices"][:]
        indptr = f_in["matrix"]["indptr"][:]
        shape = f_in["matrix"]["shape"][:]
        X = csc_matrix((data, indices, indptr), shape=shape)

        # Load features
        features = f_in["matrix"]["features"]
        gene_ids = features["id"][:]
        gene_ids_str = [g.decode() if isinstance(g, bytes) else g for g in gene_ids]

        # Build mask to remove ENSG genes
        mask = np.array([not g.startswith("ENSG") for g in gene_ids_str])

        # Filter matrix rows
        X_filtered = X[mask, :]

        # Save new filtered file
        with h5py.File(outfile, "w") as f_out:
            grp_matrix = f_out.create_group("matrix")
            Xcsc = csc_matrix(X_filtered)
            grp_matrix.create_dataset("data", data=Xcsc.data)
            grp_matrix.create_dataset("indices", data=Xcsc.indices)
            grp_matrix.create_dataset("indptr", data=Xcsc.indptr)
            grp_matrix.create_dataset("shape", data=Xcsc.shape)
            grp_matrix.create_dataset("barcodes", data=f_in["matrix"]["barcodes"][:])

            grp_features = grp_matrix.create_group("features")
            for key in features.keys():
                ds = features[key][:]
                if ds.shape[0] == mask.shape[0]:
                    # dataset is per-gene → apply mask
                    grp_features.create_dataset(key, data=ds[mask])
                else:
                    # dataset is global or metadata → copy as is
                    grp_features.create_dataset(key, data=ds)

# --- Batch processing with preserved folder structure ---
input_root = Path("/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes")
# output_root = Path("/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/filtered_cellranger_h5/Pf_all_genes")
output_root = Path("/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs_rmHsapiens/Pf_all_genes")

print("Starting ENSG gene removal from h5 files...")

# Loop over both raw and filtered feature files
for infile in input_root.glob("5736STDY*/outs/*_feature_bc_matrix.h5"):
    # Compute relative path
    rel_path = infile.relative_to(input_root)
    print(f"Processing {infile}...")

    # Create corresponding output folder
    out_dir = output_root / rel_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keep the same filename (so YASCP recognises it)
    outfile = out_dir / infile.name

    # Call your filtering function
    filter_h5_ensg(infile, outfile)
    print(f"Filtered {infile} -> {outfile}")
    
    # Copy metrics_summary.csv if it exists
    csv_file = infile.parent / "metrics_summary.csv"
    
    shutil.copy2(csv_file, out_dir / csv_file.name)
    print(f"Copied {csv_file.name} to {out_dir}")
   
