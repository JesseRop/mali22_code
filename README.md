# Mali 2022 P. falciparum single-Cell RNA-seq analysis

## Overview

This repository contains code used to analyze approximately 20 asymptomatic and symptomatic *Plasmodium falciparum* infections as part of a collaboration between the [Lawniczak lab](https://www.sanger.ac.uk/group/lawniczak-group/) at the Wellcome Sanger Institute and the Djimdes lab in the MRTC (Malaria Research and Training Center) Mali. This analysis was conducted as part of a PhD project at the University of Cambridge.

Samples were processed using 10X Genomics technology as per previously published protocols and papers. Key outcomes include stage-specific transcriptomics, strain identification and differential expression (DE) between strains, and comparative analysis between symptomatic and asymptomatic infections

## Analysis summary

This pipeline is designed to process 10X Genomics data from raw sequencing reads through cell calling, quality control, stage annotation, strain deconvolution, strain DE, and asymptomatic VS symptomatic analysis. Each numbered directory represents a sequential step in the analysis pipeline.

## Pipeline overview

### Stage 0: Initial preprocessing
**Directory**: `0_initial_preprocessing/`

Preparation of sample metadata and initial file processing, including creation of decode files for sample tracking.

**Key scripts:**
- `m22_decode_file_prep_jumpcode.Rmd` - Prepare sample tracking and metadata files
- `m2_sympto_mdata.Rmd` - Generate symptomatic/asymptomatic metadata tables
- `m2_sympto_mdata_manuscript_tables.Rmd` - Format metadata for manuscript

### Stage 1: Cell Ranger processing
**Directory**: `1_cellranger/`

Align sequencing reads to the *P. falciparum* reference genome and perform initial cell barcode calling using Cell Ranger.

**Key scripts:**
- `cellranger_count_pf_allgenes_hs_lsf_*.sh` - Submit Cell Ranger count jobs to LSF (includes human genome reference for doublet detection)
- `m2_remove_HsGns_from_raw.Rmd` - Remove human genes from raw expression matrices for downstream analysis

**Note**: Cell Ranger configuration includes both *P. falciparum* (main focus) and *Homo sapiens* reference genomes for improved technical doublet detection.

### Stage 2: Cell calling & Quality Control
**Directory**: `2_cell_calling/`

Identify true cells vs. ambient RNA using multiple cell calling approaches and quality assessment.

**Key tools/scripts:**
- **EmptyDrops**: `edrops.nf` - Nextflow pipeline for EmptyDrops analysis
- **CellBender**: `cbender_main.nf`, `cbender_main_default.nf` - Nextflow pipeline for CellBender ambient RNA removal
- **Assessment**: `cellbender_learning_curve_analysis.Rmd` - Evaluate cell calling performance
- Supporting scripts: `get_estimates_from_umi_counts.py`, `cbender_rm_human_gns_frm_cranger.py`

**Output**: Quality metrics and cell/background probability estimates

### Stage 3: Preliminary stage annotation
**Directory**: `3_prelim_stage_annotation/`

Assign preliminary parasite developmental stage annotations using reference-based annotation methods.

**Key scripts:**
- `m2_prelim_stg_annotn.Rmd` - Preliminary stage annotation pipeline
- `m2_subset_raw_2_CbendrCr.Rmd` - Subset data to cell-called cells
- `submsn_scmap.nf` / `submsn_scmap_multi_ref.nf` - scmap for stage cell type annotation
- `submsn_singleR.nf` - SingleR for automated stage annotation
- `singleR_anotn_comprsns_btwn_CbenderCr_iterations.Rmd` - Compare annotation approaches

### Stage 4: Preliminary QC after cell calling
**Directory**: `4_prelim_QC_aftr_cell_calling/`

Initial quality control assessment following cell identification.

**Key scripts:**
- `m21_raw_QC_prelim_cbendr.Rmd` - Comprehensive QC report on CellBender-processed data

### Stage 5: Strain deconvolution (Souporcell)
**Directory**: `5_souporcell/`

Identify *P. falciparum* strain composition within each infection using genomic variants.

**Key scripts:**
- `soupc_v25_m2.nf` - Main Souporcell deconvolution pipeline
- `soupc_v25_m2_k5_verif.nf` - Verification with k=5 clusters
- `soupc_v25_m2_lr.nf` - Long-read variant integration
- `m2_strain_decon_all.Rmd` - Strain deconvolution analysis and visualization

**Output**: Strain assignments for each cell

### Stage 6: Doublet QC
**Directory**: `6_doublet_QC/`

Identify and remove technical doublets (two cells captured in one droplet) and doublet-like clusters.

**Key tools/scripts:**
- **DoubletFinder**: `doubletF_rand_clust.nf` - Nextflow wrapper for DoubletFinder
- Analysis: `m21_raw_QC_prelim_cbendr_dblts.Rmd` - Doublet detection and visualization

### Stage 7: Pseudobulk genotyping & pre-QC
**Directory**: `7_pseudobulk_gtyping_preQC/`

Generate pseudobulk genotypes from single-cell data for strain verification and assessment.

**Key scripts:**
- `pbulk_gtyping.nf` - Generate pseudobulk genotypes from cell-level variants
- `m2_strain_deconv_refnd_all.Rmd` - Refined strain deconvolution analysis
- `m2_strain_deconv_K_verification_all.Rmd` - Validate strain cluster numbers
- `strain_pbulk_gtyping_supp_data_n_figs.Rmd` - Supplementary figures for strain analysis

### Stage 8: Remove strain doublets & poor QC Cells
**Directory**: `8_remove_strain_doublets_poorstrnQC/`

Filter out cells with poor strain typing quality and strain-level doublets.

**Key scripts:**
- `remove_strain_doublets_and_doublet_like_clusters.Rmd` - Final doublet and poor QC filtering

### Stage 9: Clean stage annotation & final QC
**Directory**: `9_cln_stage_annotation_final_QC/`

Refine stage annotations and apply final quality control filtering to generate analysis-ready dataset.

**Key scripts:**
- `m22_qcd_anotn.Rmd` - Clean and finalize stage annotations
- `m2_anotd_strnDblts_stgUnlabld_rm.Rmd` - Remove unlabeled and doublet stages
- `m2_low_female_investigation.Rmd` - Investigate low female gametocyte representation

### Stage 10: Integration
**Directory**: `10_integration/`

Integrate single-cell data across samples using batch-correction methods (Seurat/Harmony integration, SCT).

**Key scripts:**
- `m22_raw_integration.Rmd` - Full dataset integration
- `m22_raw_integration_asex.Rmd` - Asexual-stage specific integration
- `m2_integration_replace_BPCells_counts.Rmd` - BPCells optimization for large datasets

### Stage 11: Integrated asexual exploration
**Directory**: `11_integrated_asex_exploration/`

Exploratory analysis of integrated asexual-stage parasites across all samples.

**Key scripts:**
- `m22_raw_integration_asex_QC_clustrs.Rmd` - Cluster quality assessment
- `m22_raw_integration_asex_pseudotime.Rmd` - Trajectory inference and pseudotime analysis
- `commited_sct_clusters_trajectory.Rmd` - Committed cluster trajectory analysis
- `m2_integration_asx_clustr_gn_explrtn.Rmd` - Gene expression exploration within clusters

### Stage 12: Integrated asexual MCA & Lasse integration
**Directory**: `12_integrated_asex_mca_lasse/`

Multiple correspondence analysis (MCA) and integration with external reference datasets (Lasse datasets).

### Stage 13: Stage-specific Differential Expression (DE) & markers
**Directory**: `13_stage_DE_markers/`

Identify stage-specific marker genes and perform differential expression analysis.

**Key tools/scripts:**
- `m2_stg_markers.Rmd` - Stage marker identification
- `stg_MAST_markers_volcano_plots.Rmd` - MAST differential expression with visualization
- `GO_cnmf_gep.Rmd` - Gene set analysis and cNMF analysis
- `GO_cnmf_gep_manuscript_supp_data.Rmd` - Manuscript supplementary data generation
- `m2_scvi_to_seurat_afm.Rmd` - scVI integration and conversion to Seurat
- Supporting Python: `cnmf.py` - cNMF (consensus non-negative matrix factorization) implementation

### Stage 14: Strain-specific DE
**Directory**: `14_strain_DE/`

Compare gene expression between different strains within shared developmental stages.

### Stage 15: Asymptomatic vs. symptomatic comparison
**Directory**: `15_asymptomatic_vs_symptomatic/`

Main comparative analysis between asymptomatic and symptomatic infections.

**Key scripts:**
- `m2_sympto_DE.Rmd` - Symptomatic vs. asymptomatic differential expression
- `m2_sympto_mdata_asex.Rmd` - Metadata comparison across clinical phenotypes

### Additional Analysis modules

**GO Analysis** - `GO_analysis/`
- `GO_gaf_dbase_generator.Rmd` - Generate GO annotation databases
- `go_prep.R` - GO preparation utilities

**IBD Estimation** - `IBD_estimation/`
- Identity by descent analysis

**Host Exploration** - `host_exploration/`
- Analysis of human host cell contamination and responses

**Supporting Resources**

- `pf_common_vars.R` - Shared color schemes, stage definitions, and analysis parameters across all scripts
- `project_plotting_fns.R` - Common plotting functions
- `bash_scripts/` - Collection of helper shell scripts for pipeline submission and data processing

### Genome references
- **PlasmoDB** (*P. falciparum* Genome annotations)
- **NCBI Gene Ontology**: Functional annotation

## Project metadata

- **Author**: Jesse Rop (Lawniczak Lab, Wellcome Sanger Institute)
- **License**: MIT License (see [LICENSE](LICENSE) file)
- **Collaborators**: Djimdes Lab, MRTC Mali
- **Institution**: University of Cambridge, Wellcome Sanger Institute

## File organization

- Numbered directories (`0_*` through `15_*`): Sequential analysis stages
- R Markdown files (`.Rmd`): Documented analysis workflows with embedded code
- Nextflow scripts (`.nf`): Automated pipeline workflows
- Shell scripts (`.sh`): LSF cluster submission scripts
- Python scripts/notebooks: Specialized analyses (cNMF, etc.)
- `html_notebooks/` / `farm_rmd_nbooks/`: Rendered analysis outputs


## Citation

If you use code from this repository, please cite the associated publication and acknowledge the Lawniczak and Djimdes laboratories.

---

**Last Updated**: January 2026
