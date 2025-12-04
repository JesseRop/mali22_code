## VARIABLES

rank_colours <- c("CellEd" = "red", "Non-Cell" = "black", "Below_Limit" = "grey60")


# FUNCTIONS
## Barcode ranks plotting function incorporating thresholds for the UMI cut off for the knee and inflection for the number of likely cells and total droplets calculated by droplet utils - Added as geom hlines

bcranks_cutoff_plot_fn <-  function(dset, hline_droplets = "umi_counts_cutoff_total_droplets_", hline_cells = "umi_counts_cutoff_cell_estimate_", hline_cells_lab = "umi_cutoff_cells:", hline_droplets_lab = "umi_cutoff_droplets:", line_part = "knee", ncols = 6) {
  label_data = distinct(dset, donor, species, jumpcode, !!rlang::sym(paste0(hline_cells, line_part)), !!rlang::sym(paste0(hline_droplets, line_part)))
  
  dset %>%
    ggplot(., aes(x = rank, y = total)) +
    geom_point(size = 0.3, alpha = 0.6) +
    # scale_colour_manual(values = rank_colours) +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(aes(yintercept = !!rlang::sym(paste0(hline_cells, line_part)))) +
    geom_hline(aes(yintercept = !!rlang::sym(paste0(hline_droplets, line_part)))) +
    geom_text(data = label_data,
              aes(x = 500, 
                  y = !!rlang::sym(paste0(hline_cells, line_part)), 
                  label = paste(hline_cells_lab, !!rlang::sym(paste0(hline_cells, line_part)))),
              vjust = -0.2) +
    geom_text(data = label_data,
              aes(x = 500, 
                  y = !!rlang::sym(paste0(hline_droplets, line_part)), 
                  label = paste(hline_droplets_lab, !!rlang::sym(paste0(hline_droplets, line_part)))),
              vjust = -0.2) +
    facet_wrap2(vars(donor, species, jumpcode), ncol = ncols) +
    labs(title = paste0(line_part)) +
    theme_classic()
  
  
}



## Barcode ranks plotting function incorporating thresholds for the number of likely cells and total droplets calculated by droplet utils - Added as geom vlines

bcranks_ncells_ndrops_plot_fn <-  function(dset, don_nm, ncols = 6) {
  label_data = distinct(dset, donor, species, jumpcode, n_cells_cell_estimate_estimated_ncells, n_cells_total_droplets_estimated_ndroplets)
  
  dset %>%
    data.frame() %>%
    # filter(!(duplicated(rank))) %>%
    ggplot(., aes(x = rank, y = total)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_x_log10() +
    scale_y_log10() +
    geom_vline(aes(xintercept = n_cells_cell_estimate_estimated_ncells)) +
    geom_vline(aes(xintercept = n_cells_total_droplets_estimated_ndroplets)) +
    # labs(title = don_nm) +
    geom_text(data = label_data,
              aes(y = 5, 
                  x = n_cells_cell_estimate_estimated_ncells, 
                  label = paste0("cell_estimate: ", n_cells_cell_estimate_estimated_ncells)),
              vjust = -0.2, hjust =0.2, angle = 90) +
    geom_text(data = label_data,
              aes(y = 5, 
                  x = n_cells_total_droplets_estimated_ndroplets, 
                  label = paste0("total_droplets: ", n_cells_total_droplets_estimated_ndroplets)),
              vjust = -0.2, hjust =0.2, angle = 90) +
    facet_wrap2(vars(donor, species, jumpcode), ncol = ncols, scales = "free_x") +
    theme_classic()
  
  
  
}

