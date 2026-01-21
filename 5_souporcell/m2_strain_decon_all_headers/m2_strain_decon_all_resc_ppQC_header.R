

sample_set = 'Mali2'
# algn_nm = 'minmap'
algn_nm = 'hsat'
species = 'Pf'

soup_dir = "soupc_resc"
sv_dir = "pbulk_gtypes_ppQC_resc"

strain_levels_univrs <- c(paste0("SC", 1:20), "Doublet", "Negative")


## Declare variables

##Copy souporcell files for the optimum clusters

## Donor IDs
## All sample ids
sample_name_nms <- sort(
  c(
    "MSC33",
    "MSC48",
    "MSC49",
    "MSC55",
    "MSC60",
    "MSC66",
    # "MSC38",
    # "MSC40",
    # "MSC50_MACS",
    # "MSC53",
    # "MSC57",
    # "MSC24",
    # "MSC37",
    # "MSC39",
    # "MSC45",
    # "MSC50_SLO",
    # "MSC54",
    # "MSC67",
    "MSC68",
    # "MSC70",
    # "MSC41",
    # "MSC51",
    # "MSC25",
    # "MSC1",
    "MSC3",
    "MSC13",
    "MSC14"#,
    # "MSC12"
  )
)

## poor quality samples
sample_rm <- c("MSC1", "MSC24", "MSC37", "MSC38", "MSC39", "MSC41", "MSC51", "MSC25")

## K=1 samples
sample_k1 <- c("MSC45", "MSC50_MACS", "MSC53", "MSC54", "MSC57", "MSC67", "MSC70")

## Remove poor quality samples and those with k == 1
# sample_name_nms <- sample_name_nms[!(sample_name_nms %in% c(sample_rm, sample_k1))]
sample_name_nms <- sample_name_nms[!(sample_name_nms %in% c(sample_rm))]


## Alignment methods
algn_nms = c('minmap', 'hsat')

## Combine donor ids with alignment methods since each donor is ran once per alignment method
sample_name_nms_al <- rep(sample_name_nms, length(algn_nms))
algn_nms_al <- rep(algn_nms, each = length(sample_name_nms))

sm_aln_nms_al <- paste(rep(sample_name_nms, length(algn_nms)),
                       rep(algn_nms, each = length(sample_name_nms)),
                       sep = '_')

## number of optimum clusters to investigate per sample
# opt_clusters_nms <- list(
#   "MSC33" = c(6:10),
#   "MSC48" = c(3:7),
#   "MSC49" = c(5:9),
#   # "MSC55" = c(4:8),
#   # "MSC60" = c(5:9),
#   "MSC55" = c(5:9),
#   "MSC60" = c(6:10),
#   "MSC66" = c(7:11),
#   "MSC38" = c(4:8),
#   "MSC40" = c(6:10),
#   "MSC50_MACS" = c(2:6),
#   "MSC53" = c(1:3),
#   "MSC57" = c(1:3),
#   "MSC24" = c(1:3),
#   "MSC37" = c(5:9),
#   "MSC39" = c(2:6),
#   "MSC45" = c(1:3),
#   "MSC50_SLO" = c(2:6),
#   "MSC54" = c(1:3),
#   "MSC67" = c(1:3),
#   "MSC68" = c(2:6),
#   "MSC70" = c(1:3),
#   "MSC41" = c(2:6),
#   "MSC51" = c(1:5),
#   "MSC25" = c(2:6),
#   "MSC1" = c(3:7),
#   "MSC3" = c(1:4),
#   "MSC13" = c(1:4),
#   "MSC14" = c(5:9),
#   "MSC12" = c(1:5)
# ) %>% .[sample_name_nms]

# ## number of optimum clusters to investigate per sample
opt_clusters_nms <- list(
  "MSC33" = c(7:9),
  "MSC48" = c(3:5),
  "MSC49" = c(5:7),
  # "MSC55" = c(4:8),
  # "MSC60" = c(5:9),
  "MSC55" = c(5:8),
  "MSC60" = c(6:8),
  "MSC66" = c(7:9),
  "MSC38" = c(4:8),
  "MSC40" = c(7:9),
  "MSC50_MACS" = c(2:4),
  "MSC53" = c(1:3),
  "MSC57" = c(1:3),
  "MSC24" = c(1:3),
  "MSC37" = c(5:9),
  "MSC39" = c(2:6),
  "MSC45" = c(1:3),
  "MSC50_SLO" = c(2:4),
  "MSC54" = c(1:3),
  "MSC67" = c(1:3),
  "MSC68" = c(2:4),
  "MSC70" = c(1:3),
  "MSC41" = c(2:6),
  "MSC51" = c(1:5),
  "MSC25" = c(2:6),
  "MSC1" = c(3:7),
  "MSC3" = c(1:3),
  "MSC13" = c(1:3),
  "MSC14" = c(6:8),
  "MSC12" = c(1:3)
) %>% .[sample_name_nms]

opt_narrow_clusters_nms <-  list(
  "MSC33" = c(8:9),
  "MSC48" = c(3:5),
  "MSC49" = c(5:7),
  "MSC55" = c(6:8),
  "MSC60" = c(7:8), 
  "MSC66" = c(8:9), 
  "MSC38" = c(5:7), 
  "MSC40" = c(8:9), 
  "MSC50_MACS" = c(2), 
  "MSC53" = c(2), 
  "MSC57" = c(2), 
  "MSC24" = c(2), 
  "MSC37" = c(5:7), 
  "MSC39" = c(4:6), 
  "MSC45" = c(2), 
  "MSC50_SLO" = c(3:4), 
  "MSC54" = c(2), 
  "MSC67" = c(2), 
  "MSC68" = c(3:4), 
  "MSC70" = c(2), 
  "MSC41"= c(3:5), 
  "MSC51"= c(2), 
  "MSC25"= c(2:3), 
  "MSC1" = c(4:6), 
  "MSC3" = c(2:3), 
  "MSC13" = c(2:3), 
  "MSC14" = c(7:8), 
  "MSC12" = c(2:3)
  ) %>% .[sample_name_nms]

opt_narrow_clusters_nms_aln <- rep(opt_narrow_clusters_nms, length(algn_nms)) %>% set_names(sm_aln_nms_al)

## Select most likely optimum k
optimum_k = c("MSC33" =8, 
              "MSC48" =4, 
              "MSC49" =5, 
              "MSC55" =7, 
              "MSC60" =8, 
              "MSC66" =8, 
              "MSC38" =4, 
              "MSC40" =9, 
              "MSC50_MACS" =1, 
              "MSC53" =1, 
              "MSC57" =1, 
              "MSC24" =1, 
              "MSC37" =5, 
              "MSC39" =4, 
              "MSC45" =1, 
              "MSC50_SLO" =2, 
              "MSC54" =1, 
              "MSC67" =1, 
              "MSC68" =3, 
              "MSC70" =1, 
              "MSC41"= 3, 
              "MSC51"= 1, 
              "MSC25"= 2, 
              "MSC1" = 5, 
              "MSC3" = 2, 
              "MSC13" = 2, 
              "MSC14" = 8, 
              "MSC12" = 2) %>% .[sample_name_nms]

optimum_k_aln = rep(optimum_k, length(algn_nms)) %>% set_names(sm_aln_nms_al)


## number of optimum clusters (K) iterations for each donor
opt_clusters_len <- map(opt_clusters_nms, length)

## number of optimum clusters to investigate per sample factoring the number of alignments
opt_clusters_nms_aln <- rep(opt_clusters_nms, length(algn_nms)) %>% set_names(sm_aln_nms_al)
opt_clusters_len_aln <- map(opt_clusters_nms_aln, length)


## Donors factoring the number of alignments
sample_name <- map2(sample_name_nms, opt_clusters_len[sample_name_nms], ~
                      rep(.x, each = .y))

sample_name <- unlist(flatten(sample_name))

opt_clusters <- unlist(flatten(opt_clusters_nms[sample_name_nms]))

## sorted donor_opt_cluster names
sorted_names <- paste0(sample_name, "_", opt_clusters)

## sorted factoring alignment
sorted_names_al <- paste0(sample_name,
                          "_",
                          rep(algn_nms, each = length(sample_name)),
                          "_",
                          opt_clusters)

iter = opt_clusters_len

downsample_cells = "yes"

# ## Code for runing on farm for pseudobulk genotyping
# nextflow ../mali22_code_lnk/nf_scripts/pbulk_gtyping.nf \
# --sv_dir "pbulk_gtypes_ppQC_resc" \
# --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc_resc/hsat/parent/possorted_genome_bam.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/{MSC12,MSC13,MSC14,MSC3,MSC33,MSC40,MSC48,MSC49,MSC50_SLO,MSC55,MSC60,MSC66,MSC68}/pbulk_gtypes_ppQC_resc/bcodes/hsat/stage_afm*/*tsv" \
# --odir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" \
# --scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/cr_subset_bam_linux.sh" \
# --ref  "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"  \
# -resume

