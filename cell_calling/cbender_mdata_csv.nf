#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

params.input_h5 = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/cbender_custom_*/*/cb.h5' // change as needed
// params.input_h5 = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC3/cbender_custom_wHsPf/cb_LR0.0001_E150/cb.h5' // change as needed

// params.python = '/usr/bin/python3' // path to python; change if using env/module
params.scrpt = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/python_nbooks/cbender_mdata_csv.py'
params.obs_fname = 'cb_mdata.csv'
params.var_fname = 'cb_gene_mdata.csv'
// params.workers = 1

// Channel: find files, attach original parent dir as value
h5_ch = Channel
          .fromPath(params.input_h5, checkIfExists: true)
          .map { file -> tuple(file, file.getParent().toString(), params.obs_fname, params.var_fname) }

h5_ch.view()

// Process: run converter for each file
process CBENDER_CONVERT {
  
  tag { h5.getName() }
  cpus 1
  memory '20 GB'
  time '1.h'

  // publish outputs to the input file's parent directory (dynamic)
  // mode can be 'copy', 'symlink', or 'move'
  publishDir { orig_parent }, mode: 'copy', pattern: "*.csv", overwrite: true

  input:
    tuple path(h5), val(orig_parent), val(obs_fname), val(var_fname)

  // declare files that the process will produce in its work dir
  output:
    tuple path("${obs_fname}"), path("${var_fname}")

  script:
  """
  set -euo pipefail

  # Activate the cellbender GPU virtualenv for this job
  source /software/team222/jr35/cellbender_gpu/cellbender_gpu_env/bin/activate
  
  python ${params.scrpt} --input "${h5}" --outdir . --obs-fname "${obs_fname}" --var-fname "${var_fname}"
  """
}


workflow {
  h5_csv = CBENDER_CONVERT(h5_ch)
  h5_csv.view()
}

