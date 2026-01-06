mrm_processed_data <- sample_spectrums

# --- Shiny app for compound review ----

library(shiny)
library(tidyverse)
library(ggplot2)
library(htmltools)

# --- DATA LOADING PLACEHOLDER ---
# Ensure 'mrm_processed_data' is in your environment before running.
# mrm_processed_data <- sample_spectrums

# --- UI ---
ui <- fluidPage(
  titlePanel("Compound/Transition Review Tool"),

  # Add JavaScript for keyboard navigation
  tags$script(HTML(
    "
    $(document).on('keydown', function(e) {
      if (e.keyCode == 37) { // Left arrow
        Shiny.setInputValue('key_left', Math.random());
      } else if (e.keyCode == 39) { // Right arrow
        Shiny.setInputValue('key_right', Math.random());
      }
    });
  "
  )),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Evaluation File"),
      textInput(
        "eval_filename",
        "File Name:",
        value = "compound_evaluations.csv",
        width = "100%"
      ),
      fluidRow(
        column(
          6,
          actionButton(
            "load_csv",
            "Load Evaluations",
            class = "btn-primary",
            width = "100%"
          )
        ),
        column(
          6,
          actionButton(
            "save_csv",
            "Save Evaluations",
            class = "btn-success",
            width = "100%"
          )
        )
      ),
      hr(),
      h4("Navigation"),
      p("Use ← → arrow keys to navigate"),
      uiOutput("compound_selector"),
      uiOutput("transition_selector"), # Transition checkboxes

      hr(),
      h4("Add Evaluation"),
      textAreaInput(
        "evaluation_text",
        "Evaluation (auto-saved on navigation):",
        value = "",
        rows = 4,
        width = "100%"
      ),

      hr(),
      h4("Export Plot"),
      downloadButton(
        "download_plot",
        "Save Current Plot (.png)",
        class = "btn-info",
        width = "100%"
      ),
      downloadButton(
        "download_all_plots",
        "Save All Plots (.pdf)",
        class = "btn-warning",
        width = "100%"
      ),

      hr(),
      h5("Saved Evaluations"),
      tableOutput("saved_evaluations_table")
    ),

    mainPanel(
      width = 9,
      plotOutput("compound_plot", height = "600px")
    )
  )
)

# --- SERVER ---
server <- function(input, output, session) {
  # 1. Check Data Availability
  if (!exists("mrm_processed_data")) {
    stop(
      "Error: 'mrm_processed_data' not found in environment. Please load data."
    )
  }

  # 2. Determine Navigation Mode (Compound vs Transition)
  has_compound_names <- "compound_name" %in%
    colnames(mrm_processed_data) &&
    any(!is.na(mrm_processed_data$compound_name))

  if (has_compound_names) {
    message("Compound names found. Using compound-based navigation.")

    # Get base list
    compounds_df_base <- mrm_processed_data %>%
      filter(!is.na(compound_name)) %>%
      distinct(compound_name)

    # Get metadata from mrm_info if it exists
    if (exists("mrm_info") && is.data.frame(get("mrm_info"))) {
      mrm_info_cols <- colnames(mrm_info)
      cols_to_use <- "compound_name"

      id_col <- NA_character_
      if ("inchikey" %in% mrm_info_cols) {
        cols_to_use <- c(cols_to_use, "inchikey")
        id_col <- "inchikey"
      } else if ("cas" %in% mrm_info_cols) {
        cols_to_use <- c(cols_to_use, "cas")
        id_col <- "cas"
      }

      mass_col <- NA_character_
      if ("exact_mass" %in% mrm_info_cols) {
        cols_to_use <- c(cols_to_use, "exact_mass")
        mass_col <- "exact_mass"
      } else if ("precursor_mz" %in% mrm_info_cols) {
        cols_to_use <- c(cols_to_use, "precursor_mz")
        mass_col <- "precursor_mz"
      }

      mrm_details <- mrm_info %>%
        select(all_of(unique(cols_to_use))) %>%
        group_by(compound_name) %>%
        summarise(
          inchikey = if (!is.na(id_col)) {
            dplyr::first(!!sym(id_col))
          } else {
            NA_character_
          },
          exact_mass = if (!is.na(mass_col)) {
            dplyr::first(!!sym(mass_col))
          } else {
            NA_real_
          }
        ) %>%
        ungroup() %>%
        distinct()

      compounds_df <- left_join(
        compounds_df_base,
        mrm_details,
        by = "compound_name"
      )
    } else {
      compounds_df <- compounds_df_base %>%
        mutate(inchikey = NA_character_, exact_mass = NA_real_)
    }

    # Ensure columns exist even if join failed
    if (!"inchikey" %in% colnames(compounds_df)) {
      compounds_df$inchikey <- NA_character_
    }
    if (!"exact_mass" %in% colnames(compounds_df)) {
      compounds_df$exact_mass <- NA_real_
    }

    compounds_df <- compounds_df %>%
      arrange(compound_name) %>%
      mutate(display_name = compound_name)
  } else {
    message("No compound names found. Using transition-based navigation.")
    compounds_df <- mrm_processed_data %>%
      filter(!is.na(transition)) %>%
      distinct(transition, precursor_mz) %>%
      arrange(transition) %>%
      mutate(
        compound_name = transition,
        display_name = transition,
        inchikey = NA_character_,
        exact_mass = precursor_mz
      )
  }

  # 3. Reactive Values for Evaluations
  evaluations_data <- reactiveVal({
    tibble(
      compound_name = character(),
      inchikey = character(),
      evaluation = character(),
      timestamp = character()
    )
  })

  # Helper function to save
  save_evaluation_if_new <- function(
    compound_name,
    inchikey_val,
    evaluation_text
  ) {
    trimmed_eval <- trimws(evaluation_text)
    if (nchar(trimmed_eval) == 0) {
      return()
    }

    current_data <- evaluations_data()
    evaluation_exists <- any(
      current_data$compound_name == compound_name &
        current_data$evaluation == trimmed_eval
    )

    if (!evaluation_exists) {
      new_evaluation <- tibble(
        compound_name = compound_name,
        inchikey = inchikey_val,
        evaluation = trimmed_eval,
        timestamp = as.character(Sys.time())
      )
      evaluations_data(bind_rows(current_data, new_evaluation))
    }
  }

  # 4. Inputs
  output$compound_selector <- renderUI({
    label <- if (has_compound_names) {
      "Select Compound:"
    } else {
      "Select Transition:"
    }
    selectInput(
      "selected_compound",
      label,
      choices = setNames(compounds_df$compound_name, compounds_df$display_name),
      selected = compounds_df$compound_name[1]
    )
  })

  output$transition_selector <- renderUI({
    req(input$selected_compound, has_compound_names)
    available_transitions <- mrm_processed_data %>%
      filter(compound_name == input$selected_compound, !is.na(transition)) %>%
      distinct(transition) %>%
      pull(transition) %>%
      sort()

    if (length(available_transitions) > 1) {
      checkboxGroupInput(
        "selected_transitions",
        "Filter Transitions:",
        choices = available_transitions,
        selected = available_transitions
      )
    } else {
      NULL
    }
  })

  # 5. Data Retrieval for Current Selection
  current_compound_data <- reactive({
    req(input$selected_compound)
    if (has_compound_names) {
      mrm_processed_data %>% filter(compound_name == input$selected_compound)
    } else {
      mrm_processed_data %>% filter(transition == input$selected_compound)
    }
  })

  current_compound_info <- reactive({
    req(input$selected_compound)
    compounds_df %>% filter(compound_name == input$selected_compound)
  })

  # 6. Auto-save on Navigation
  previous_compound <- reactiveVal(NULL)

  observe({
    req(input$selected_compound)
    prev_comp <- previous_compound()

    # Save previous if changed
    if (!is.null(prev_comp) && prev_comp != input$selected_compound) {
      prev_comp_info <- compounds_df %>%
        filter(compound_name == prev_comp) %>%
        slice(1)
      save_evaluation_if_new(
        prev_comp,
        prev_comp_info$inchikey,
        input$evaluation_text
      )
    }
    previous_compound(input$selected_compound)

    # Load existing for current
    current_evals <- evaluations_data() %>%
      filter(compound_name == input$selected_compound) %>%
      arrange(desc(timestamp))

    updateTextAreaInput(
      session,
      "evaluation_text",
      value = if (nrow(current_evals) > 0) current_evals$evaluation[1] else ""
    )
  })

  # 7. Plot Logic
  current_plot_object <- reactive({
    compound_data <- current_compound_data()
    compound_info <- current_compound_info()

    if (nrow(compound_data) == 0) {
      return(ggplot() + labs(title = "No data") + theme_minimal())
    }

    all_transitions <- compound_data %>%
      filter(!is.na(transition)) %>%
      distinct(transition) %>%
      pull(transition)
    n_total_transitions <- length(all_transitions)

    # Filter Logic
    filtered_plot_data <- compound_data
    if (has_compound_names && n_total_transitions > 1) {
      if (!is.null(input$selected_transitions)) {
        filtered_plot_data <- compound_data %>%
          filter(transition %in% input$selected_transitions)
      }
    }

    if (nrow(filtered_plot_data) == 0) {
      return(
        ggplot() + labs(title = "No transitions selected") + theme_minimal()
      )
    }

    n_transitions_plotted <- length(unique(filtered_plot_data$transition[
      !is.na(filtered_plot_data$transition)
    ]))

    # Subtitles
    if (has_compound_names) {
      base_subtitle <- if (n_total_transitions > 1) {
        paste0(
          "Showing ",
          n_transitions_plotted,
          " of ",
          n_total_transitions,
          " transitions"
        )
      } else {
        paste("Number of transitions:", n_total_transitions)
      }

      # --- UPDATED SUBTITLE LOGIC ---
      # Check if mrm_info exists AND has 'exact_mass' column
      has_exact_mass_source <- exists("mrm_info") &&
        "exact_mass" %in% colnames(mrm_info)

      val <- compound_info$exact_mass[1]

      # Only show value if source was exact_mass AND value exists
      # Otherwise show NA
      mass_string <- if (has_exact_mass_source && !is.na(val)) {
        paste0(" (Exact Mass: ", sprintf("%.4f", val), ")")
      } else {
        " (Exact Mass: NA)"
      }

      final_subtitle <- paste0(base_subtitle, mass_string)
      plot_title <- paste0("Compound: ", input$selected_compound)
    } else {
      # Fallback mode (Transition navigation)
      final_subtitle <- paste0(
        "Precursor m/z: ",
        sprintf("%.4f", compound_info$exact_mass[1])
      )
      plot_title <- paste0("Transition: ", input$selected_compound)
    }

    # Build Plot
    p <- filtered_plot_data %>%
      filter(!is.na(transition)) %>%
      ggplot(aes(rtime, intensity, color = sample_name)) +
      geom_line(alpha = 0.7, linewidth = 0.5) +
      scale_color_viridis_d(option = "C") +
      labs(
        title = plot_title,
        subtitle = final_subtitle,
        x = "Retention Time (s)",
        y = "Intensity",
        color = "Sample"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold")
      ) +
      guides(
        color = guide_legend(override.aes = list(alpha = 1, linewidth = 2))
      )

    if (has_compound_names && n_transitions_plotted > 1) {
      p <- p + facet_wrap(~transition, ncol = 1, scales = "free_y")
    }
    p
  })

  output$compound_plot <- renderPlot({
    current_plot_object()
  })

  # 8. Keyboard Logic
  observeEvent(input$key_left, {
    req(input$selected_compound)
    idx <- which(compounds_df$compound_name == input$selected_compound)
    if (idx > 1) {
      updateSelectInput(
        session,
        "selected_compound",
        selected = compounds_df$compound_name[idx - 1]
      )
    }
  })

  observeEvent(input$key_right, {
    req(input$selected_compound)
    idx <- which(compounds_df$compound_name == input$selected_compound)
    if (idx < nrow(compounds_df)) {
      updateSelectInput(
        session,
        "selected_compound",
        selected = compounds_df$compound_name[idx + 1]
      )
    }
  })

  # 9. File Handlers
  observeEvent(input$load_csv, {
    req(input$eval_filename)
    if (file.exists(input$eval_filename)) {
      tryCatch(
        {
          raw <- read_csv(
            input$eval_filename,
            show_col_types = FALSE,
            col_types = cols(.default = col_character())
          )
          # Renaming fallbacks
          if ("comment" %in% names(raw) && !"evaluation" %in% names(raw)) {
            raw <- rename(raw, evaluation = comment)
          }
          if ("cas" %in% names(raw) && !"inchikey" %in% names(raw)) {
            raw <- rename(raw, inchikey = cas)
          }

          # Ensure cols exist
          req_cols <- c("compound_name", "inchikey", "evaluation", "timestamp")
          for (c in req_cols) {
            if (!c %in% names(raw)) raw[[c]] <- NA_character_
          }

          evaluations_data(raw %>% select(all_of(req_cols)))
          showNotification(
            paste("Loaded", nrow(raw), "evaluations."),
            type = "message"
          )

          # Update current text
          curr <- raw %>%
            filter(compound_name == input$selected_compound) %>%
            arrange(desc(timestamp))
          updateTextAreaInput(
            session,
            "evaluation_text",
            value = if (nrow(curr) > 0) curr$evaluation[1] else ""
          )
        },
        error = function(e) {
          showNotification(paste("Error loading:", e$message), type = "error")
        }
      )
    } else {
      showNotification("File not found.", type = "error")
    }
  })

  observeEvent(input$save_csv, {
    file_to_save <- input$eval_filename
    req(file_to_save)

    # Save current text box first
    if (!is.null(input$selected_compound) && input$evaluation_text != "") {
      info <- current_compound_info()
      save_evaluation_if_new(
        input$selected_compound,
        info$inchikey[1],
        input$evaluation_text
      )
    }

    tryCatch(
      {
        write_csv(evaluations_data(), file_to_save)
        showNotification(paste("Saved to", file_to_save), type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving:", e$message), type = "error")
      }
    )
  })

  # 10. Downloads
  output$download_plot <- downloadHandler(
    filename = function() {
      safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", input$selected_compound)
      paste0("plot_", safe_name, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(
        file,
        plot = current_plot_object(),
        device = "png",
        width = 11,
        height = 8,
        dpi = 300
      )
    }
  )

  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("all_compound_plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      showNotification(
        "Generating PDF... please wait.",
        id = "plot_prog",
        duration = NULL
      )
      pdf(file, width = 11, height = 8)

      for (i in seq_len(nrow(compounds_df))) {
        cmp <- compounds_df$compound_name[i]
        cmp_info <- compounds_df[i, ]

        # Get Data
        p_data <- if (has_compound_names) {
          mrm_processed_data %>% filter(compound_name == cmp)
        } else {
          mrm_processed_data %>% filter(transition == cmp)
        }
        if (nrow(p_data) == 0) {
          next
        }

        # Recreate Plot Logic (simplified for PDF loop)
        n_trans <- length(unique(p_data$transition[!is.na(p_data$transition)]))

        if (has_compound_names) {
          subt <- paste0("Total transitions: ", n_trans)

          # --- UPDATED PDF LOGIC ---
          has_exact_mass_source <- exists("mrm_info") &&
            "exact_mass" %in% colnames(mrm_info)
          mass_val <- cmp_info$exact_mass

          if (has_exact_mass_source && !is.na(mass_val)) {
            mass_str <- sprintf("%.4f", mass_val)
          } else {
            mass_str <- "NA"
          }

          subt <- paste0(subt, " (Exact Mass: ", mass_str, ")")
          ti <- paste0("Compound: ", cmp)
        } else {
          subt <- paste0("Precursor: ", sprintf("%.4f", cmp_info$exact_mass))
          ti <- paste0("Transition: ", cmp)
        }

        p <- p_data %>%
          filter(!is.na(transition)) %>%
          ggplot(aes(rtime, intensity, color = sample_name)) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          scale_color_viridis_d(option = "C") +
          labs(
            title = ti,
            subtitle = subt,
            x = "Retention Time (s)",
            y = "Intensity"
          ) +
          theme_bw() +
          theme(legend.position = "bottom")

        if (has_compound_names && n_trans > 1) {
          p <- p + facet_wrap(~transition, ncol = 1, scales = "free_y")
        }
        print(p)
      }
      dev.off()
      removeNotification("plot_prog")
      showNotification("PDF Saved!", type = "message")
    }
  )

  output$saved_evaluations_table <- renderTable(
    {
      req(input$selected_compound)
      evaluations_data() %>%
        filter(compound_name == input$selected_compound) %>%
        arrange(desc(timestamp)) %>%
        select(timestamp, evaluation) %>%
        slice(1:3)
    },
    width = "100%"
  )
}

# --- Run ---
shinyApp(ui, server)
