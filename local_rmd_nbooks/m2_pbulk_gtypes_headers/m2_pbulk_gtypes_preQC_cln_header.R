
## Sample names with more than one strain
# sample_name_nms <- c("MSC33","MSC37","MSC38","MSC39","MSC40", "MSC48", "MSC49", "MSC55", "MSC60", "MSC66", "MSC68", "MSC50_SLO", "MSC25", "MSC41") %>% sort()
sample_name_nms <- sort(c("MSC33", "MSC48","MSC49", "MSC55", "MSC60", "MSC66", "MSC38", "MSC40", "MSC50_MACS", "MSC53", "MSC57", "MSC24", "MSC37", "MSC39", "MSC45", "MSC50_SLO", "MSC54", "MSC67", "MSC68", "MSC70", "MSC41", "MSC51", "MSC25", "MSC1", "MSC3", "MSC13", "MSC14"))

## Remove MSC38 & MSC39
sample_name_nms <- sample_name_nms[!(sample_name_nms %in% c("MSC37","MSC38","MSC39","MSC1"))]
sample_name_nms <- sort(c("MSC55", "MSC60"))

# Get the optimum QCd K
opt_narrow_clusters_nms <-  list(
  "MSC33" = c(8),
  "MSC37" = c(5),
  "MSC38" = c(4),
  "MSC39" = c(4),
  "MSC40" = c(8),
  "MSC48" = c(4),
  "MSC49" = c(6),
  "MSC55" = c(6),
  "MSC60" = c(7),
  "MSC66" = c(8),
  "MSC68" = c(3),
  "MSC50_SLO" = c(2),
  "MSC50_MACS" = c(1),
  "MSC53" = c(1),
  "MSC57" = c(1),
  "MSC24" = c(1),
  "MSC45" = c(1),
  "MSC54" = c(1),
  "MSC67" = c(1),
  "MSC70" = c(1),
  "MSC41" = c(3),
  "MSC51" = c(1),
  "MSC25" = c(2),
  "MSC1" = c(5),
  "MSC3" = c(2),
  "MSC13" = c(2),
  "MSC14" = c(8)
  
) %>% .[sample_name_nms]


sample_name_k <- map2(sample_name_nms,
                      opt_narrow_clusters_nms,
                      ~ rep(.x, each = length(.y))) %>% unlist()

# Path donor label
sample_name_nms_dir <- sample_name_nms

# Path k label
opt_narrow_clusters_nms_dir <-  opt_narrow_clusters_nms

## sorted donor_opt_cluster names
sorted_names <- imap(opt_narrow_clusters_nms, ~ paste0(.y, "_", .x)) %>% unlist(., recursive =
                                                                                  F, use.names = F)

optimum_k <- opt_narrow_clusters_nms %>% unlist() %>% set_names(sorted_names)

# Get the pseudobulk groups/categories
grps <- c('strain', 'stage_ag_strain')

##Alignment incorporation
# algn_nms = c('minmap', 'hsat')
algn_nms = c('minmap')
algn_nm = 'minmap'

# Directory to read pseudobulk genotype results from
gtyp_step = "pbulk_gtypes_preQC_cln"

sample_set <- "Mali2"

species <- "Pf"
mt_chrom_nm <- "Pf3D7_MIT_v3"

## NB-DONT DELETE BELOW- RUN ONCE Get genome size like this cut -f1,2 fasta.fai > genome.size
# system("cut -f1,2 /lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta.fai > /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pdb66_Genome.size")

chrm_sz <- "../../Pf3D7_genomes_gtfs_indexs/Pdb66_Genome.size"
chrm_prfx <- "Pf3D7_0|Pf3D7_|_v3"

## Directory to save barcodes for the next iteration of pseudobulk genotyping
sv_dir = "pbulk_gtypes_preQC_cln_final"

# nextflow nf_scripts/pbulk_gtyping.nf --sv_dir "pbulk_gtypes_preQC_cln" --save_bam_sset true --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc/*/parent/possorted_genome_bam.bam" --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/bcodes/minmap/{stage_ag_strain*,strain*}/*tsv" --odir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" --mapper "minmap" --scrpt "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh" --ref  "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"