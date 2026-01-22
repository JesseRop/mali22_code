#!/usr/bin/env Rscript
# Extract CellBender learning curve assessment from HTML reports
#
# This script parses CellBender cb_report.html files to extract the automated
# learning curve assessment, which indicates whether the learning rate was
# appropriate or needs adjustment.
#
# Author: jr35
# Date: 24 November 2025

library(tidyverse)
library(rvest)
library(glue)


#' Extract learning curve assessment from CellBender HTML report
#'
#' @param html_file Path to the cb_report.html file
#' @return Character string with the assessment, or NA if not found
extract_learning_curve_assessment <- function(html_file) {
  tryCatch({
    # Read HTML
    html_content <- read_html(html_file)

    # Find all paragraphs
    paragraphs <- html_content %>% 
      html_nodes("p") %>% 
      html_text(trim = TRUE)

    # Find "Summary:" paragraph
    summary_idx <- which(paragraphs == "Summary:")

    if (length(summary_idx) > 0) {
      # Get the next few paragraphs after "Summary:"
      # The assessment should be in one of the next paragraphs
      for (i in (summary_idx + 1):min(summary_idx + 5, length(paragraphs))) {
        text <- paragraphs[i]
        if (grepl("This learning curve looks normal", text) || 
            grepl("This is unusual behavior", text) ||
            grepl("This is slightly unusual behavior", text) ||
            grepl("reduced --learning-rate", text)) {
          return(text)
        }
      }
    }

    # Alternative: search all paragraphs for the patterns
    for (text in paragraphs) {
      if (grepl("This learning curve looks normal", text) || 
          grepl("This is unusual behavior", text) ||
          grepl("This is slightly unusual behavior", text) ||
          grepl("reduced --learning-rate", text)) {
        return(text)
      }
    }

    return(NA_character_)

  }, error = function(e) {
    message(glue("Error processing {html_file}: {e$message}"))
    return(NA_character_)
  })
}

#' Extract and save the learning curve image from CellBender HTML report
#'
#' @param html_file Path to the cb_report.html file
#' @param out_png Path to save the PNG file
#' @return TRUE if image saved, FALSE otherwise
extract_learning_curve_image <- function(html_file, out_png) {
  tryCatch({
    html_content <- read_html(html_file)
    # Target the first image AFTER the convergence heading
    target_img <- html_content %>% html_node(xpath = "//h2[@id='Assessing-convergence-of-the-algorithm']/following::img[1]")
    img_src <- if (!is.null(target_img)) html_attr(target_img, "src") else NA_character_
    # Fallback: previous heuristic (first base64 image) if heading-based extraction fails
    if (is.na(img_src) || !grepl("^data:image/png;base64,", img_src)) {
      img_nodes <- html_content %>% html_nodes("img")
      img_srcs <- img_nodes %>% html_attr("src")
      img_idx <- which(grepl("^data:image/png;base64,", img_srcs))
      if (length(img_idx) == 0) return(FALSE)
      img_src <- img_srcs[img_idx[1]]
    }
    img_b64 <- sub("^data:image/png;base64,", "", img_src)
    img_bin <- base64enc::base64decode(img_b64)
    writeBin(img_bin, out_png)
    return(TRUE)
  }, error = function(e) {
    message(glue("Error extracting image from {html_file}: {e$message}"))
    return(FALSE)
  })
}


#' Parse file path to extract sample name and parameters
#'
#' @param file_path Full path to the cb_report.html file
#' @return Named list with sample, learning_rate, and epochs
parse_file_path <- function(file_path) {
  # Extract sample name
  sample <- str_extract(file_path, "MSC[^/]+")
  if (is.na(sample)) sample <- "Unknown"
  
  # Extract learning rate
  lr_match <- str_extract(file_path, "cb_LR([0-9.]+)")
  learning_rate <- str_replace(lr_match, "cb_LR", "")
  if (is.na(learning_rate)) learning_rate <- "Unknown"
  
  # Extract epochs
  epochs_match <- str_extract(file_path, "_E(\\d+)")
  epochs <- str_replace(epochs_match, "_E", "")
  if (is.na(epochs)) epochs <- "Unknown"
  
  list(
    sample = sample,
    learning_rate = learning_rate,
    epochs = epochs
  )
}


#' Categorize the assessment text
#'
#' @param assessment The assessment text from the HTML
#' @return Category string
categorize_assessment <- function(assessment) {
  if (is.na(assessment)) {
    return("Not Found")
  } else if (grepl("This learning curve looks normal", assessment)) {
    return("Normal")
  } else if (grepl("This is unusual behavior, and a reduced --learning-rate is indicated", assessment)) {
    return("Needs Lower LR")
  } else if (grepl("This is slightly unusual behavior, and a reduced --learning-rate might be indicated", assessment)) {
    return("Maybe Needs Lower LR")
  } else if (grepl("unusual behavior", assessment) || grepl("reduced --learning-rate", assessment)) {
    return("Maybe Needs Lower LR")
  } else {
    return("Other")
  }
}


# Main script
main <- function(search_pattern = NULL, output_dir = NULL) {
  # Define the search pattern
  if (is.null(search_pattern)) {
    # Default pattern
    base_path <- "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf"
    pattern <- file.path(base_path, "MSC*/cbender_custom2_rmHs/cb_LR*/cb_report.html")
  } else {
    pattern <- search_pattern
  }
  
  # Define output directory
  if (is.null(output_dir)) {
    # Default output directory
    output_dir <- "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/cell_calling"
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    message(glue("Creating output directory: {output_dir}"))
    dir.create(output_dir, recursive = TRUE)
  }
  
  message(glue("Searching for files matching: {pattern}"))
  
  # Find all matching files
  html_files <- Sys.glob(pattern)
  
  if (length(html_files) == 0) {
    stop(glue("No files found matching pattern: {pattern}"))
  }
  
  message(glue("Found {length(html_files)} HTML files to process\n"))
  
  # Extract condition from pattern for output filename (move up so it's in scope)
  condition <- if (is.null(search_pattern)) {
    "rmHs"
  } else if (grepl("cbender_custom2_rmHs", search_pattern)) {
    "rmHs"
  } else if (grepl("cbender_custom2_wHsPf", search_pattern)) {
    "wHsPf"
  } else {
    # Extract from pattern or use generic
    gsub(".*/cbender_custom2_([^/]+)/.*", "\\1", search_pattern, perl = TRUE)
  }

  # Process each file
  results <- map_dfr(html_files, function(html_file) {
    message(glue("Processing: {html_file}"))

    # Extract path information
    path_info <- parse_file_path(html_file)

    # Extract assessment
    assessment <- extract_learning_curve_assessment(html_file)

    # Categorize
    category <- categorize_assessment(assessment)

    # Save learning curve image (PNG) for each run, with condition in filename
    img_outfile <- file.path(output_dir, glue("{path_info$sample}_LR{path_info$learning_rate}_E{path_info$epochs}_{condition}_learning_curve.png"))
    img_saved <- extract_learning_curve_image(html_file, img_outfile)
    if (img_saved) {
      message(glue("  Saved learning curve image: {img_outfile}"))
    } else {
      message(glue("  No learning curve image found in: {html_file}"))
    }

    # Return as data frame row
    tibble(
      Sample = path_info$sample,
      Learning_Rate = path_info$learning_rate,
      Epochs = path_info$epochs,
      Assessment_Category = category,
      Full_Assessment = ifelse(is.na(assessment), "N/A", assessment),
      File_Path = html_file,
      Learning_Curve_Image = ifelse(img_saved, img_outfile, NA_character_)
    )
  })
  
  # Sort results
  results <- results %>%
    arrange(Sample, Learning_Rate, Epochs)
  
  # ...existing code...
  
  # Save to CSV
  output_file <- file.path(output_dir, glue("cellbender_learning_curve_assessments_{condition}.csv"))
  write_csv(results, output_file)
  
  message("\n", strrep("=", 80))
  message(glue("Results saved to: {output_file}"))
  message(strrep("=", 80), "\n")
  
  # Print summary statistics
  message("Summary Statistics:")
  message(glue("Total files processed: {nrow(results)}\n"))
  
  message("Assessment Categories:")
  print(table(results$Assessment_Category))
  
  message("\nBreakdown by Sample:")
  sample_summary <- results %>%
    count(Sample, Assessment_Category) %>%
    pivot_wider(names_from = Assessment_Category, values_from = n, values_fill = 0)
  
  print(sample_summary)
  
  # Save summary table
  summary_file <- file.path(output_dir, glue("cellbender_learning_curve_summary_{condition}.csv"))
  write_csv(sample_summary, summary_file)
  message(glue("\nSummary table saved to: {summary_file}"))
  
  return(results)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required packages
if (!requireNamespace("rvest", quietly = TRUE)) {
  message("Installing required package: rvest")
  install.packages("rvest")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  message("Installing required package: tidyverse")
  install.packages("tidyverse")
}
if (!requireNamespace("base64enc", quietly = TRUE)) {
  message("Installing required package: base64enc")
  install.packages("base64enc")
}

# Run the main function with optional arguments
if (length(args) > 0) {
  # User provided a pattern
  search_pattern <- args[1]
  output_dir <- if (length(args) > 1) args[2] else NULL
  
  message(glue("\nUsing user-provided pattern: {search_pattern}"))
  if (!is.null(output_dir)) {
    message(glue("Using user-provided output directory: {output_dir}"))
  }
  message("")
  
  results <- main(search_pattern, output_dir)
} else {
  # Use default pattern and output directory
  message("\nNo pattern provided, using default (cbender_custom2_rmHs)\n")
  message("Usage: Rscript extract_cellbender_learning_curve_assessment.R [pattern] [output_dir]\n")
  message("Example patterns:")
  message('  Rscript extract_cellbender_learning_curve_assessment.R "/path/to/MSC*/cbender_custom2_rmHs/cb_LR*/cb_report.html"')
  message('  Rscript extract_cellbender_learning_curve_assessment.R "/path/to/MSC*/cbender_custom2_wHsPf/cb_LR*/cb_report.html" "/path/to/output"')
  message("")
  
  results <- main()
}
