# Contrast example ----

## packages and functions ----

source("r/lrms_analysis_scripts.R")

## fns ----

# identifies mzML files in input directory, imports and makes a df with spectral info and metadata for each sample
# since
process_srm_batch <- function(input, mrm_info = NULL, mz_tolerance = 0.15) {
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

    # Standardize column names
    colnames(gc_df)[4] <- "precursor_mz"
    colnames(gc_df)[8] <- "product_mz"
    colnames(gc_df)[7] <- "ce"

    # Clean names
    gc_df <- gc_df %>% clean_names()

    # PERFORM FUZZY JOIN AT THE FEATURE LEVEL (much faster!)
    if (!is.null(mrm_info)) {
      # Ensure mrm_info has the required columns
      required_cols <- c("compound_name", "precursor_mz", "product_mz")
      if (!all(required_cols %in% colnames(mrm_info))) {
        stop(
          "mrm_info must contain columns: compound_name, precursor_mz, product_mz"
        )
      }

      # Perform fuzzy join on small feature dataframe (67 rows instead of millions)
      # Perform fuzzy join on small feature dataframe (67 rows instead of millions)
      gc_df_annotated <- gc_df %>%
        difference_inner_join(
          mrm_info %>%
            rename(
              precursor_mz_ref = precursor_mz,
              product_mz_ref = product_mz
            ),
          by = c(
            "product_mz" = "product_mz_ref",
            "precursor_mz" = "precursor_mz_ref"
          ),
          max_dist = mz_tolerance
        ) %>%
        # Add transition label
        mutate(
          transition = paste0(
            round(precursor_mz, digits = 1),
            " -> ",
            round(product_mz, digits = 1)
          )
        )
    } else {
      # If no annotation, add NA columns
      gc_df_annotated <- gc_df %>%
        mutate(
          compound_name = NA_character_,
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

    # Now extract chromatogram data only for annotated transitions
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
  ) %>%
    mutate(
      tidy_data = map(
        path,
        ~ {
          cat("Processing:", basename(.x), "\n")
          process_single_srm(.x)
        }
      )
    ) %>%
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

  blank_colors <- rep("black", n_blanks) # All blanks are black
  standard_colors <- rep("red", n_standards) # All standards are red

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

## Processing ----

# mrms for annotating transitions
mrm_info <- read_csv(
  "/Users/a33241/Documents/GitHub/reproducible_ms_lowres/files/251103_mrms_rev5.csv"
)

# creating dataframe with spectral and chromatographic info by sample
tic()
sample_spectrums <- process_srm_batch(
  input = "/Users/a33241/HRMS_data/contrast/GC/251103_rev4mrms/",
  # mrm_info = NULL
  read_csv("files/251103_mrms_rev5.csv")
) %>%
  # Append sample metadata from file name
  mutate(
    # Extract temperature (number before c_)
    temperature = as.numeric(str_extract(file_name, "\\d+(?=c_)")),
    # Extract ce (number before ce)
    ce = as.numeric(str_extract(file_name, "\\d+(?=ce)")),
    # Create sample_name combining temperature and energy
    sample_name = paste0(temperature, "C_", ce, "eV"),
    # Convert to ordered factor by temperature first, then ce
    sample_name = factor(
      sample_name,
      levels = paste0(
        rep(c(200, 285), each = 4), # each temperature repeated 4 times
        "C_",
        rep(c(5, 23, 42, 60), times = 2), # all 4 ce levels per temperature
        "eV"
      ),
      ordered = TRUE
    )
  )
toc()


sample_spectrums %>% glimpse
