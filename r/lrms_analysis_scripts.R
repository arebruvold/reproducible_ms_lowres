# Load the necessary libraries ----
library(MSnbase)
library(tidyverse)
library(janitor)
library(fuzzyjoin)
library(tictoc)


ms_metadataer <- function(ms_data) {
  # Get list of mzML files in the directory
  ms_files <- list.files(ms_data, pattern = "\\.mzML$", full.names = TRUE)

  # Check if any files were found
  if (length(ms_files) == 0) {
    stop("No mzML files found in directory: ", ms_data)
  }

  # Open the first mzML file (suppress warnings from openMSfile)
  ms_file <- suppressWarnings(openMSfile(ms_files[[1]]))

  # Get instrument and run info (suppress warnings from these functions)
  instrument <- suppressWarnings(instrumentInfo(ms_file))
  run <- suppressWarnings(runInfo(ms_file))

  # Return as a list instead of printing
  list(
    instrument = instrument,
    run = run
  )
}


process_srm_batch <- function(input, mrm_info = NULL, mz_tolerance = 0.1) {
  # Load RColorBrewer for color palette
  library(RColorBrewer)

  # Determine if input is a directory or file(s)
  if (length(input) == 1 && dir.exists(input)) {
    # Input is a directory, get all mzML files
    file_paths <- list.files(
      path = input,
      pattern = "(?i)\\.mzml$",
      full.names = TRUE
    )

    if (length(file_paths) == 0) {
      stop("No mzML files found in directory: ", input)
    }

    cat("Found", length(file_paths), "mzML files in directory:", input, "\n")
  } else if (all(file.exists(input))) {
    # Input is file path(s)
    file_paths <- input
    cat("Processing", length(file_paths), "specified file(s)\n")
  } else {
    stop("Input must be either a valid directory path or existing file path(s)")
  }

  # Helper function to process a single SRM file
  process_single_srm <- function(file_path) {
    # Import raw data
    gc_raw <- readSRMData(file_path)

    # Feature info
    gc_df <- fData(gc_raw)

    # Add SRM index column
    gc_df$srm_index <- row(gc_df)[, 1]

    # Standardize column names FIRST, before any other operations
    colnames(gc_df)[4] <- "precursor_mz"
    colnames(gc_df)[8] <- "product_mz"
    colnames(gc_df)[7] <- "ce"

    # PERFORM FUZZY JOIN AT THE FEATURE LEVEL (much faster!)
    if (!is.null(mrm_info)) {
      # Ensure mrm_info has the required columns
      required_cols <- c("compound_name", "precursor_mz", "product_mz")
      if (!all(required_cols %in% colnames(mrm_info))) {
        stop(
          "mrm_info must contain columns: compound_name, precursor_mz, product_mz"
        )
      }

      # Perform fuzzy join on small feature dataframe
      gc_df_annotated <- gc_df %>%
        difference_inner_join(
          mrm_info %>%
            dplyr::rename(
              precursor_mz_ref = precursor_mz,
              product_mz_ref = product_mz
            ),
          by = c(
            "product_mz" = "product_mz_ref",
            "precursor_mz" = "precursor_mz_ref"
          ),
          max_dist = mz_tolerance
        ) %>%
        # Keep the original precursor_mz and product_mz from gc_df (with .x suffix)
        mutate(
          # precursor_mz = precursor_mz.x,
          # product_mz = product_mz.x,
          transition = paste0(
            round(precursor_mz, digits = 1),
            " -> ",
            round(product_mz, digits = 1)
          )
        ) %>%
        select(
          -ends_with(".x"),
          -ends_with(".y"),
          precursor_mz,
          product_mz,
          everything()
        )
    } else {
      # If no annotation, add NA columns and keep all transitions
      gc_df_annotated <- gc_df %>%
        mutate(
          compound_name = NA_character_,
          precursor_mz_ref = NA_real_,
          product_mz_ref = NA_real_,
          cas = NA_character_,
          transition = paste0(
            round(precursor_mz, digits = 1),
            " -> ",
            round(product_mz, digits = 1)
          )
        )
      cat(
        "No compound annotation provided; proceeding without annotation.\n"
      )
    }

    # Now extract chromatogram data for all transitions (annotated or not)
    gc_tib <- gc_df_annotated %>%
      mutate(
        chrom = map(
          srm_index,
          ~ tibble(
            rtime = rtime(gc_raw[.x]) * 60, # Convert to seconds
            intensity = intensity(gc_raw[.x])
          )
        )
      ) %>%
      unnest(chrom) %>%
      select(-matches("offset")) # Remove offset columns

    return(gc_tib)
  }

  # Process all files
  cat("Processing", length(file_paths), "mzML files...\n")

  all_data <- tibble(
    path = file_paths,
    file_name = basename(file_paths)
  ) |>
    mutate(
      tidy_data = map(
        path,
        \(.x) {
          cat("Processing:", basename(.x), "\n")
          process_single_srm(.x)
        }
      )
    ) |>
    unnest(tidy_data)

  # ADD SAMPLE TYPE CLASSIFICATION BASED ON FILE NAMES
  # Define regex patterns for sample types
  blank_regex <- "(?i)(blank|blk|rinse|wash)"
  standard_regex <- "(?i)(std|standard|cal|calibr)"

  # Classify sample types based on file names
  all_data <- all_data %>%
    mutate(
      sample_type = case_when(
        str_detect(file_name, blank_regex) ~ "blank",
        str_detect(file_name, standard_regex) ~ "standard",
        TRUE ~ "sample"
      )
    )

  # CREATE COLOR PALETTE BASED ON SAMPLE TYPES
  # Count samples by type
  sample_counts <- all_data %>%
    distinct(file_name, sample_type) %>%
    count(sample_type)

  n_samples <- sample_counts$n[sample_counts$sample_type == "sample"]
  n_blanks <- sample_counts$n[sample_counts$sample_type == "blank"]
  n_standards <- sample_counts$n[sample_counts$sample_type == "standard"]

  # Handle cases where some sample types might not exist
  if (length(n_samples) == 0) {
    n_samples <- 0
  }
  if (length(n_blanks) == 0) {
    n_blanks <- 0
  }
  if (length(n_standards) == 0) {
    n_standards <- 0
  }

  # Create color vectors
  sample_colors <- if (n_samples > 0) {
    RColorBrewer::brewer.pal(max(n_samples, 3), "Set1")[1:n_samples]
  } else {
    character(0)
  }

  blank_colors <- rep("black", n_blanks)
  standard_colors <- rep("red", n_standards)

  # Create a lookup table for file colors
  file_colors <- all_data %>%
    distinct(file_name, sample_type) %>%
    arrange(sample_type, file_name) %>%
    mutate(
      row_idx = row_number(),
      color = case_when(
        sample_type == "sample" ~ sample_colors[pmin(
          row_idx,
          length(sample_colors)
        )],
        sample_type == "blank" ~ "black",
        sample_type == "standard" ~ "red"
      )
    ) %>%
    select(file_name, sample_type, color)

  # Add colors to the main dataset
  all_data <- all_data %>%
    left_join(file_colors, by = c("file_name", "sample_type"))

  cat("Sample type classification:\n")
  cat("- Samples:", n_samples, "\n")
  cat("- Blanks:", n_blanks, "\n")
  cat("- Standards:", n_standards, "\n")

  if (!is.null(mrm_info)) {
    cat(
      "Successfully annotated",
      length(unique(all_data$compound_name[!is.na(all_data$compound_name)])),
      "unique compounds\n"
    )
  }

  return(all_data)
}
