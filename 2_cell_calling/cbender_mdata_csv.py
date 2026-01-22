#!/usr/bin/env python3
"""
cbender_mdata_csv.py

Simple, robust converter that reads a Cellbender `cb.h5` (anndata) and writes:
  - obs -> cb_mdata.csv (cell metadata)
  - var -> cb_gene_mdata.csv (feature/gene metadata)

Designed to be called from a wrapper script (e.g. cbender_submsn_cellNo_edropsNo.sh)
right after Cellbender completes. Example usage in the wrapper:

If --outdir is omitted or set to 'auto'|'same'|'.' outputs go to the input file's folder.
Call from wrapper, e.g.:
  python cbender_mdata_csv.py --input "$o_file"
  python cbender_mdata_csv.py --input "$o_file" --outdir "$o_dir"
  
Requirements: Cellbender environment with
  cellbender.remove_background.downstream.anndata_from_h5
  source /software/team222/jr35/cellbender_gpu/cellbender_gpu_env/bin/activate

"""
from __future__ import annotations
import argparse
import glob
import os
import sys
from pathlib import Path
import pandas as pd

# Try to import the downstream converter from CellBender (standard env)
try:
    from cellbender.remove_background.downstream import anndata_from_h5
except Exception:
    anndata_from_h5 = None


def parse_args():
    p = argparse.ArgumentParser(description="Convert Cellbender cb.h5 anndata to CSV (obs & var).")
    p.add_argument("--input", required=True, help="File or glob to .h5 (e.g. '/path/*/cb.h5').")
    p.add_argument("--outdir", required=False, default=None,
                   help="Output directory or template; may include {sample}. If omitted or 'auto' outputs go to input folder.")
    p.add_argument("--obs-fname", default="cb_mdata.csv")
    p.add_argument("--var-fname", default="cb_gene_mdata.csv")
    p.add_argument("--sample", help="Optional sample name for templating {sample} in --outdir.")
    p.add_argument("--follow-symlinks", action="store_true", default=True)
    p.add_argument("--overwrite", action="store_true", default=False)
    return p.parse_args()


def write_from_h5(h5_path: str, outdir: str, obs_fname: str, var_fname: str, follow_symlinks: bool, overwrite: bool):
    if anndata_from_h5 is None:
        raise RuntimeError("cellbender downstream not available. Activate correct environment.")
    real = os.path.realpath(h5_path) if follow_symlinks else h5_path
    if not os.path.exists(real):
        raise FileNotFoundError(real)
    adata = anndata_from_h5(real)
    outp = Path(outdir)
    outp.mkdir(parents=True, exist_ok=True)
    obs_path = outp / obs_fname
    var_path = outp / var_fname
    if not overwrite and obs_path.exists() and var_path.exists():
        print(f"Skipping (exists): {real} -> {outp}")
        return
    adata.obs.to_csv(obs_path)
    adata.var.to_csv(var_path)
    print(f"Wrote: {obs_path}, {var_path}")


def main():
    args = parse_args()
    files = sorted(glob.glob(args.input))
    if not files:
        print(f"No files match: {args.input}", file=sys.stderr); return 2
    for f in files:
        sample = args.sample or Path(f).parent.name or Path(f).stem
        # Determine outdir: explicit -> use (allow {sample}), 'auto'|'same'|'.' or omitted -> input parent
        if args.outdir is None or str(args.outdir).lower() in ("auto", "same", "."):
            outdir = str(Path(f).parent)
        else:
            outdir = args.outdir.format(sample=sample) if "{sample}" in args.outdir else args.outdir
        # if outdir appears to be a file path, use its parent dir
        if outdir.endswith(".csv"):
            outdir = str(Path(outdir).parent)
        try:
            write_from_h5(f, outdir, args.obs_fname, args.var_fname, args.follow_symlinks, args.overwrite)
        except Exception as e:
            print(f"ERROR processing {f}: {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

