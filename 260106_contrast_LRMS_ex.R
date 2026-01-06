# Contrast example ----

## packages and functions ----

source("r/lrms_analysis_scripts.R")


## Processing ----

# mrms for annotating transitions
mrm_info <- read_csv(
  "/Users/a33241/Documents/GitHub/reproducible_ms_lowres/files/251103_mrms_rev5.csv"
)

# directory with MS data files
ms_data <- "ms_data/"

# get metadata/info

ms_metadataer(ms_data)

# creating dataframe with spectral and chromatographic info by sample
{
  tic()
  sample_spectrums <- process_srm_batch(
    input = ms_data,
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
}


## Visualization ----

sample_spectrums %>%
  filter(compound_name == "Bis(2-ethylhexyl)hexanedioate") %>%
  ggplot(aes(rtime, intensity, color = sample_name)) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  scale_color_viridis_d(option = "C") +
  labs(
    title = "Bis(2-ethylhexyl)hexanedioate",
    x = "Retention Time (seconds)",
    y = "Intensity",
    color = "Sample"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, linewidth = 2))
  ) +
  facet_wrap(~transition)
