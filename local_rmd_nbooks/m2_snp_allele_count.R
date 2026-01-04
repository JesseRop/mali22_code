# R
# For HbS and HbC
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)
library(conflicted)

conflict_prefer("filter", "dplyr")

## MSC14
# nextflow  run get_snp_alleles.nf --bam "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY11771546/outs/possorted_genome_bam.bam" --postn "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code/nf_scripts/pfsa_pos.txt"

## farm nextflow script for samtools mpile up
# nextflow run get_snp_alleles.nf --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/strain_k*/new_bams/*bam"
# nextflow run get_snp_alleles.nf --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/strain_k*_minmap/new_bams/*bam" --postn "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code/nf_scripts/pfsa_pos.txt" -resume
# nextflow run get_snp_alleles.nf --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_preQC_cln/strain_k*_minmap/new_bams/*bam" --postn "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code/nf_scripts/who_art_pos.txt" -resume


# Set the directory containing the text files
# mydir <- "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"
# file_list <- list.files(path = mydir, pattern = "snp\\.txt$", full.names = TRUE, recursive = T) %>%
#   set_names(nm = (str_remove_all(.,'/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf//|/snp.txt') %>% str_replace_all(., c( '/'='_')) )) 
mydir_14 <- "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC14"

file_list_14 <- list.files(path = mydir_14, pattern = "snp\\.txt$", full.names = TRUE, recursive = T) %>%
  set_names(nm = (str_remove_all(.,'/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/|/snp_don/snp.txt') %>% str_replace_all(., c( '/'='_')) )) 

mydir <- "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf"

file_list_m2 <- list.files(path = mydir, pattern = ".*_snp\\.txt$", full.names = TRUE, recursive = T) %>%
  set_names(nm = (str_remove_all(.,'/lustre.*new_bams_snp/|_snp.txt') %>% str_replace_all(., c( '/'='_')) ))  

file_list <- c(file_list_14, file_list_m2)

file_list_pfsa_who <- list.files(path = mydir, pattern = ".*_snp\\.txt$", full.names = TRUE, recursive = T) %>%
  set_names(nm = (str_remove_all(.,'/lustre.*new_bams_snp/|_snp.txt') %>% str_replace_all(., c( '/'='_')) ))  %>%
  # .[str_detect(., "pfsa|who")] %>%
  .[str_detect(., "pfsa")]

file_list_pfsa_who <- c(file_list_14, file_list_pfsa_who)

snps <- file_list_pfsa_who %>%
  map(., possibly(~read.table(.x, header = FALSE, sep = "\t")), "no supporting reads")


# file_names<- list.files(path = mydir, pattern = "*.txt", full.names = FALSE)
# file_names <- gsub("\\.txt$", "",file_names)
# my_data <- list()
# for (i in seq_along(file_list)) {my_data[[i]] <- read.table(file_list[i], header = FALSE, sep = "\t") }
my_data <- do.call(rbind, snps)




myset <- c("a","t","c","g","A","T","C","G")
df <- data.frame(matrix(nrow = nrow(my_data), ncol = 8))
colnames(df) <- myset
for (i in 1:length(myset)) { df[,myset[i]] <- str_count(my_data$V5, myset[i]) }

my_data <- cbind(my_data[, c(paste0("V", 1:4))],df)

my_data_long <- my_data %>% rownames_to_column("donor") %>% mutate(across(donor, ~str_remove(.,"\\.[0-9]$"))) %>% pivot_longer(cols = c("a", "t", "c","g",  "A", "T", "C", "G"), names_to = "nucleotide", values_to = "count") %>% filter(count >0) 
# rnames <- rep(names(file_list), each = 2)
# rnames <- ifelse(sequence(length(rnames)) %% 2 == 0, paste0(rnames, "_2"), rnames)
# rownames(my_data) <- rnames

# write.csv(my_data,"frequency_HBB_pileups_HBB_SC.csv")

my_data_long %>% group_by(donor, V2) %>% mutate(prop = count/sum(count)) %>% ggplot(aes(x = donor, y =prop, fill = nucleotide)) + geom_col() + facet_wrap(V2 ~ .) + coord_flip()
my_data_long %>% ggplot(aes(x = donor, y =count, fill = nucleotide)) + geom_col() + facet_wrap(V2 ~ .) + coord_flip()

my_data_long %>% ggplot(aes(x = donor, y =count, fill = nucleotide)) + geom_col() + facet_wrap(V2 ~ ., nrow = 1, scales = "free_x") + coord_flip()

