#!/usr/bin/env python3
"""
Filter out human genes (ENSG*) from CellRanger outputs:
- h5 matrices (*_feature_bc_matrix.h5)
- MTX matrices (raw_feature_bc_matrix, filtered_feature_bc_matrix)
"""

from pathlib import Path
import argparse
import h5py
import numpy as np
import pandas as pd
import gzip
import shutil
import csv
from scipy.sparse import csc_matrix
from scipy.io import mmread, mmwrite

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
        
        # Count genes
        total_genes = len(gene_ids_str)
        retained_genes = mask.sum()
        removed_genes = total_genes - retained_genes

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
    
    return total_genes, retained_genes, removed_genes

def filter_mtx_ensg(src_dir, dest_dir, filter_prefix="ENSG"):
    """Filter ENSG genes from MTX format (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)"""
    required = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
    if not all((src_dir / f).exists() for f in required):
        return None, None, None
    
    # Read and filter features
    features = pd.read_csv(src_dir / "features.tsv.gz", sep="\t", header=None)
    features.columns = ["id", "name", "type"]
    mask = ~features["id"].str.startswith(filter_prefix)
    kept_features = features[mask].reset_index(drop=True)
    
    # Filter matrix
    mtx = mmread(src_dir / "matrix.mtx.gz").tocsr()
    filtered_mtx = mtx[mask.values, :]
    
    # Write filtered outputs
    dest_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest_dir / "features.tsv.gz", "wt") as f:
        kept_features.to_csv(f, sep="\t", header=False, index=False)
    with gzip.open(dest_dir / "barcodes.tsv.gz", "wb") as f_out, \
         gzip.open(src_dir / "barcodes.tsv.gz", "rb") as f_in:
        f_out.write(f_in.read())
    mmwrite(dest_dir / "matrix.mtx", filtered_mtx)
    import os
    os.system(f"gzip -f {dest_dir / 'matrix.mtx'}")
    
    return len(features), len(kept_features), len(features) - len(kept_features)

def main():
    parser = argparse.ArgumentParser(description="Remove ENSG (human) genes from CellRanger h5 and MTX matrices.")
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
        help="Root folder where filtered files will be written (default: rmHsapiens Pf_all_genes)",
    )
    parser.add_argument(
        "--skip-mtx",
        action="store_true",
        help="Skip MTX format matrices (only process h5 files)",
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

    print(f"Starting ENSG gene removal for irods_id in CSV: {csv_path}")
    print(f"Found {len(irods_ids)} unique irods_id entries")

    for irods_id in irods_ids:
        # Process h5 files
        pattern = f"{irods_id}/outs/*_feature_bc_matrix.h5"
        matched = False
        for infile in input_root.glob(pattern):
            matched = True
            rel_path = infile.relative_to(input_root)
            print(f"\n[h5] Processing {infile.name} from {irods_id}...")

            out_dir = output_root / rel_path.parent
            out_dir.mkdir(parents=True, exist_ok=True)
            outfile = out_dir / infile.name

            total, retained, removed = filter_h5_ensg(infile, outfile)
            print(f"  Genes: {total} total | {retained} retained | {removed} removed (ENSG)")
            print(f"  Output: {outfile}")

            # Copy metrics_summary.csv if it exists
            csv_file = infile.parent / "metrics_summary.csv"
            if csv_file.exists():
                shutil.copy2(csv_file, out_dir / csv_file.name)
                print(f"  Copied: {csv_file.name}")
        
        if not matched:
            print(f"Warning: No h5 found for irods_id '{irods_id}'")
        
        # Process MTX files (default behavior, unless --skip-mtx flag is set)
        if not args.skip_mtx:
            sample_outs = input_root / irods_id / "outs"
            if sample_outs.exists():
                for folder in ["raw_feature_bc_matrix", "filtered_feature_bc_matrix"]:
                    src_dir = sample_outs / folder
                    if not src_dir.exists():
                        continue
                    
                    dest_dir = output_root / irods_id / "outs" / folder
                    print(f"\n[MTX] Processing {irods_id}/{folder}...")
                    
                    total, retained, removed = filter_mtx_ensg(src_dir, dest_dir)
                    if total:
                        print(f"  Genes: {total} total | {retained} retained | {removed} removed (ENSG)")
                        print(f"  Output: {dest_dir}")
                    else:
                        print(f"  ⚠️  Missing required files, skipped")

if __name__ == "__main__":
    main()
   
