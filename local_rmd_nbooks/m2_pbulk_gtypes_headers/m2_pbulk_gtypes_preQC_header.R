

## Sample names with more than one strain
# sample_name_nms <- c("MSC50_SLO", "MSC41", "MSC51", "MSC25") %>% str_sort(numeric = T)
sample_name_nms <- c("MSC33", "MSC48","MSC49", "MSC55", "MSC60", "MSC66", "MSC38", "MSC40", "MSC50_MACS", "MSC53", "MSC57", "MSC24", "MSC37", "MSC39", "MSC45", "MSC50_SLO", "MSC54", "MSC67", "MSC68", "MSC70", "MSC41", "MSC51", "MSC25") %>% sort()


# Get the optimum QCd K
opt_narrow_clusters_nms <-  list(
  "MSC33" = c(7:9),
  "MSC48" = c(3:5),
  "MSC49" = c(5:7),
  "MSC55" = c(5:7),
  "MSC60" = c(6:8),
  "MSC66" = c(7:9),
  "MSC38" = c(5:7),
  "MSC40" = c(7:9),
  "MSC50_MACS" = c(2:4),
  "MSC53" = c(2),
  "MSC57" = c(2),
  "MSC24" = c(2),
  "MSC37" = c(5:7),
  "MSC39" = c(4:6),
  "MSC45" = c(2),
  "MSC50_SLO" = c(2:4),
  "MSC54" = c(2),
  "MSC67" = c(2),
  "MSC68" = c(2:4),
  "MSC70" = c(2),
  "MSC41" = c(3:5),
  "MSC51" = c(2),
  "MSC25" = c(2:3)
) %>% .[sample_name_nms]
# optimum_k = opt_narrow_clusters_nms

sample_name_k <- map2(sample_name_nms,
                      opt_narrow_clusters_nms,
                      ~ rep(.x, each = length(.y))) %>% unlist()

## sorted donor_opt_cluster names
sorted_names <- imap(opt_narrow_clusters_nms, ~ paste0(.y, "_", .x)) %>% unlist(., recursive =
                                                                                  F, use.names = F)

optimum_k <- opt_narrow_clusters_nms %>% unlist() %>% set_names(sorted_names)

# Get the pseudobulk groups/categories
grps <- c('stage_afm_strain')

##Alignment incorporation
algn_nms = c('minmap')
algn_nm = 'minmap'

# algn_nms = c('hsat')
# algn_nm = 'hsat'

gtyp_step = "pbulk_gtypes_preQC"
# gtyp_step = "pbulk_gtypes_GE_postQC"

sample_set <- "Mali2"

species <- "Pf"
mt_chrom_nm <- "Pf3D7_MIT_v3"

## NB-DONT DELETE BELOW- RUN ONCE Get genome size like this cut -f1,2 fasta.fai > genome.size
# system("cut -f1,2 /lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta.fai > /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pdb66_Genome.size")

chrm_sz <- "../../Pf3D7_genomes_gtfs_indexs/Pdb66_Genome.size"
chrm_prfx <- "Pf3D7_0|Pf3D7_|_v3"

## Directory to save barcodes for the next iteration of pseudobulk genotyping
sv_dir = "pbulk_gtypes_preQC_cln"

## Select between downsampled and full cell barcodes. Downsampled allows faster runtime while maintaining power - yes to downsample, any other value not to downsample
downsample_cells = "yes"

## Nextflow script for pseudobulk gtyping
## Code for runing on farm for pseudobulk genotyping
##### MINIMAP
# nextflow pbulk_gtyping.nf \
# --sv_dir "pbulk_gtypes_preQC_cln" \
# --save_bam_sset true \
# --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc/*/parent/possorted_genome_bam.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/bcodes/minmap/{stage_ag_strain*,strain*}/*tsv" \
# --odir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" \
# --mapper "minmap" \
# --scrpt "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh" \
# --ref  "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

##### HISAT
# nextflow pbulk_gtyping.nf \
# --sv_dir "pbulk_gtypes_preQC_cln" \
# --save_bam_sset true \
# --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc/*/parent/possorted_genome_bam.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/bcodes/minmap/{stage_ag_strain*,strain*}/*tsv" \
# --odir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" \
# --mapper "minmap" \
# --scrpt "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh" \
# --ref  "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"



