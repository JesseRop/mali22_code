# This script reads the sample names from a CSV file and then iterates over each sample name to read the corresponding h5 file and save its metadata as a CSV file.
# Run as below
        ## source /software/team222/jr35/cellbender_gpu/cellbender_gpu_env/bin/activate
        ## python cbender_rm_human_gns_in_h5.py

# Load required modules and functions
from pathlib import Path
import argparse
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

# --- Batch processing using irods_id from CSV ---
import csv

def main():
    parser = argparse.ArgumentParser(description="Remove ENSG (human) genes from CellRanger h5 matrices based on a list of irods_id in a CSV.")
    parser.add_argument(
        "--csv",
        dest="csv_path",
        default="/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_jcode_merge_decode.csv",
        help="Path to CSV file with an 'irods_id' column (default: current project CSV)",
    )
    parser.add_argument(
        "--input-root",
        dest="input_root",
        default="/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes",
        help="Root folder containing Cell Ranger outputs (default: Pf_all_genes run root)",
    )
    parser.add_argument(
        "--output-root",
        dest="output_root",
        default="/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs_rmHsapiens/Pf_all_genes",
        help="Root folder where filtered h5 files will be written (default: rmHsapiens Pf_all_genes)",
    )
    args = parser.parse_args()

    input_root = Path(args.input_root)
    output_root = Path(args.output_root)

    csv_path = Path(args.csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    # Read irods_id values from CSV
    irods_ids = set()
    with open(csv_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        if "irods_id" not in reader.fieldnames:
            raise ValueError(f"CSV {csv_path} is missing required 'irods_id' column. Columns found: {reader.fieldnames}")
        for row in reader:
            val = (row.get("irods_id") or "").strip()
            if val:
                irods_ids.add(val)

    if not input_root.exists():
        raise FileNotFoundError(f"Input root not found: {input_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"Starting ENSG gene removal from h5 files for irods_id in CSV: {csv_path}")
    print(f"Found {len(irods_ids)} unique irods_id entries")

    for irods_id in irods_ids:
        # Some irods_id values are single, some are underscore-separated pairs
        # We expect directory names matching the irods_id exactly under input_root
        pattern = f"{irods_id}/outs/*_feature_bc_matrix.h5"
        matched = False
        for infile in input_root.glob(pattern):
            matched = True
            rel_path = infile.relative_to(input_root)
            print(f"Processing {infile}...")

            out_dir = output_root / rel_path.parent
            out_dir.mkdir(parents=True, exist_ok=True)
            outfile = out_dir / infile.name

            filter_h5_ensg(infile, outfile)
            print(f"Filtered {infile} -> {outfile}")

            # Copy metrics_summary.csv if it exists
            csv_file = infile.parent / "metrics_summary.csv"
            if csv_file.exists():
                shutil.copy2(csv_file, out_dir / csv_file.name)
                print(f"Copied {csv_file.name} to {out_dir}")
        if not matched:
            print(f"Warning: No input h5 found for irods_id '{irods_id}' using pattern '{pattern}'")

if __name__ == "__main__":
    main()
   
