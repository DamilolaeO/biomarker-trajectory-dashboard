###############################################################################
# REUSABLE BIOMARKER TRAJECTORY DASHBOARD TEMPLATE
###############################################################################
#
# SUMMARY
# This Shiny application presents aggregate results from two longitudinal
# biomarker trajectory analyses. It displays trajectory patterns, 
# aggregated class characteristics, landmark progression-free and 
# overall survival plots, and Cox regression estimates.
#
# REQUIRED FOLDER STRUCTURE
# dashboard/
#   app.R
#   dashboard_data/   # RDS files; not publicly supplied
#   www/              # static images
###############################################################################


############################################################
# 1. PACKAGES
############################################################


library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(ggplot2)
library(ggtext)


############################################################
# 2. TEMPLATE SETTINGS — EDIT THIS SECTION ONLY
############################################################

# Replace every value written in CAPITALS. The class labels must match the
# labels stored in your aggregate RDS files. Internal object names
# such as `pp` and `b2m` are retained only to keep the template code stable.


dashboard_config <- list(
  app_title = "YOUR STUDY: BIOMARKER TRAJECTORIES",
  study_name = "YOUR STUDY NAME",
  trajectory_window = "YOUR TRAJECTORY WINDOW",
  landmark_label = "YOUR LANDMARK TIME",
  adjusted_covariates = "YOUR ADJUSTMENT VARIABLES",
  stage_label = "YOUR STAGING VARIABLE",
  pathway_label = "YOUR TREATMENT GROUP VARIABLE",
  secretory_type_label = "YOUR DISEASE SUBTYPE VARIABLE",
  light_chain_label = "YOUR ADDITIONAL BIOMARKER VARIABLE",
  
  biomarker_a = list(
    name = "BIOMARKER A FULL NAME",
    short_name = "BIOMARKER A",
    cohort_n = "N = XXXX",
    class_labels = c(
      "BIOMARKER A CLASS 1",
      "BIOMARKER A CLASS 2",
      "BIOMARKER A CLASS 3",
      "BIOMARKER A CLASS 4"
    ),
    term_prefix = "trajectory_class"
  ),
  
  biomarker_b = list(
    name = "BIOMARKER B FULL NAME",
    short_name = "BIOMARKER B",
    cohort_n = "N = XXXX",
    class_labels = c(
      "BIOMARKER B CLASS 1",
      "BIOMARKER B CLASS 2",
      "BIOMARKER B CLASS 3"
    ),
    term_prefix = "trajectory_b2m_class"
  ),
  
  files = list(
    trajectory_plots = "trajectory_plots.rds",
    biomarker_a_characteristics = "biomarker_a_characteristics.rds",
    biomarker_b_characteristics = "biomarker_b_characteristics.rds",
    dashboard_outputs = "dashboard_outputs.rds",
    full_cohort_age = "full_cohort_age.png",
    full_cohort_sex = "full_cohort_sex.png",
    full_cohort_stage = "full_cohort_stage.png",
    full_cohort_pathway = "full_cohort_pathway.png",
    biomarker_a_pfs = "biomarker_a_pfs.png",
    biomarker_a_os = "biomarker_a_os.png",
    biomarker_b_pfs = "biomarker_b_pfs.png",
    biomarker_b_os = "biomarker_b_os.png"
  )
)

############################################################
# 3. LOAD PREPARED AGGREGATE OUTPUTS
############################################################

# The aggregate RDS files are intentionally not included in the public
# repository. This check gives a clear message if the app is opened without
# dashboard inputs.
dashboard_files <- c(
  trajectory_plots = file.path(
    "dashboard_data",
    dashboard_config$files$trajectory_plots
  ),
  pp_characteristics = file.path(
    "dashboard_data",
    dashboard_config$files$biomarker_a_characteristics
  ),
  b2m_characteristics = file.path(
    "dashboard_data",
    dashboard_config$files$biomarker_b_characteristics
  ),
  outputs = file.path(
    "dashboard_data",
    dashboard_config$files$dashboard_outputs
  )
)

missing_dashboard_files <- dashboard_files[
  !file.exists(dashboard_files)
]

if (length(missing_dashboard_files) > 0) {
  stop(
    paste0(
      "Approved aggregate dashboard files are not available: ",
      paste(basename(missing_dashboard_files), collapse = ", "),
      ". Patient-level data must not be added to this repository."
    )
  )
}

trajectory_plots <- readRDS(
  dashboard_files[["trajectory_plots"]]
)

pp_characteristics_summary <- readRDS(
  dashboard_files[["pp_characteristics"]]
)

b2m_characteristics_summary <- readRDS(
  dashboard_files[["b2m_characteristics"]]
)

dashboard_outputs <- readRDS(
  dashboard_files[["outputs"]]
)



############################################################
# 4. THEME
############################################################

dashboard_theme <- bs_theme(
  version = 5,
  primary = "#1565C0",
  secondary = "#607D8B",
  success = "#00897B",
  bg = "#F4F7FB",
  fg = "#17212B"
)


############################################################
# 5. TRACE NAMES
############################################################

pp_trace_names <- rep(
  dashboard_config$biomarker_a$class_labels,
  times = 2
)

b2m_trace_names <- rep(
  dashboard_config$biomarker_b$class_labels,
  times = 2
)


############################################################
# 6. HELPER FUNCTIONS
############################################################

format_p_value <- function(p_value) {
  
  if (is.na(p_value)) {
    return(NA_character_)
  }
  
  if (p_value < 0.001) {
    return("<0.001")
  }
  
  sprintf(
    "%.3f",
    p_value
  )
}


prepare_cox_display <- function(data) {
  
  data %>%
    mutate(
      term = gsub(
        dashboard_config$biomarker_a$term_prefix,
        "",
        term
      ),
      
      term = gsub(
        dashboard_config$biomarker_b$term_prefix,
        "",
        term
      ),
      
      `Hazard ratio` = sprintf(
        "%.2f",
        estimate
      ),
      
      `95% CI` = sprintf(
        "%.2f–%.2f",
        conf.low,
        conf.high
      ),
      
      `p-value` = vapply(
        p.value,
        format_p_value,
        character(1)
      )
    ) %>%
    select(
      Outcome = outcome,
      Model = model,
      `Trajectory class` = term,
      `Hazard ratio`,
      `95% CI`,
      `p-value`
    )
}

# PATIENT-CHARACTERISTICS GRAPHS

dashboard_blues <- c(
  "#0B2E59",
  "#1565C0",
  "#56A8E0",
  "#9ED2F2",
  "#B0B7C3"
)


characteristics_theme <- theme_minimal(
  base_size = 15
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 20,
      colour = "#17212B"
    ),
    
    plot.subtitle = element_text(
      size = 12,
      colour = "#607D8B",
      margin = margin(b = 18)
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 11,
      colour = "#17212B"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "#17212B"
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_blank(),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    plot.margin = margin(
      20,
      30,
      20,
      20
    )
  )

make_characteristics_plot <- function(
    summary_data,
    class_column,
    class_levels,
    characteristic
) {
  
  # AGE
  
  if (characteristic == "age") {
    
    plot_data <- summary_data$age %>%
      mutate(
        class_display = factor(
          .data[[class_column]],
          levels = class_levels
        )
      ) %>%
      filter(
        !is.na(class_display)
      )
    
    return(
      ggplot(
        plot_data,
        aes(
          x = class_display,
          fill = class_display
        )
      ) +
        geom_boxplot(
          aes(
            ymin = minimum,
            lower = lower,
            middle = middle,
            upper = upper,
            ymax = maximum
          ),
          stat = "identity",
          width = 0.6,
          alpha = 0.9
        ) +
        scale_fill_manual(
          values = dashboard_blues
        ) +
        labs(
          title = "Age at baseline",
          subtitle = paste(
            "Median and interquartile range within each",
            "trajectory class"
          ),
          x = NULL,
          y = "Age, years"
        ) +
        guides(
          fill = "none"
        ) +
        characteristics_theme
    )
  }
  
  
  # SELECT AGGREGATED CATEGORICAL DATA
  
  plot_data <- summary_data[[characteristic]]
  
  
  graph_title <- switch(
    characteristic,
    
    sex = "Sex",
    r2iss = dashboard_config$stage_label,
    pathway = dashboard_config$pathway_label,
    secretory = dashboard_config$secretory_type_label,
    light_chain = dashboard_config$light_chain_label
  )
  
  
  legend_title <- switch(
    characteristic,
    
    sex = "Sex",
    r2iss = dashboard_config$stage_label,
    pathway = dashboard_config$pathway_label,
    secretory = dashboard_config$secretory_type_label,
    light_chain = dashboard_config$light_chain_label
  )
  
  
  category_levels <- switch(
    characteristic,
    
    sex = c(
      "Female",
      "Male",
      "Missing"
    ),
    
    r2iss = c(
      "I",
      "II",
      "III",
      "IV",
      "Missing"
    ),
    
    pathway = c(
      "Transplant Eligible",
      "Transplant Not Eligible",
      "Missing"
    ),
    
    secretory = c(
      "Intact immunoglobulin",
      "Light-chain-only",
      "Non-secretor",
      "Missing"
    ),
    
    light_chain = c(
      "Kappa",
      "Lambda",
      "Not recorded",
      "Missing"
    )
  )
  
  
  plot_data <- plot_data %>%
    mutate(
      class_display = factor(
        .data[[class_column]],
        levels = class_levels
      ),
      
      category_display = factor(
        category_display,
        levels = category_levels
      ),
      
      percentage_label = if_else(
        percentage >= 7,
        paste0(
          round(percentage, 1),
          "%"
        ),
        ""
      )
    ) %>%
    filter(
      !is.na(class_display),
      !is.na(category_display)
    )
  
  
  ggplot(
    plot_data,
    aes(
      x = class_display,
      y = percentage,
      fill = category_display
    )
  ) +
    geom_col(
      width = 0.62,
      colour = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        label = percentage_label
      ),
      position = position_stack(
        vjust = 0.5
      ),
      colour = "white",
      fontface = "bold",
      size = 4.5
    ) +
    scale_fill_manual(
      values = setNames(
        dashboard_blues[
          seq_along(category_levels)
        ],
        category_levels
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = scales::label_percent(
        scale = 1
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    labs(
      title = graph_title,
      subtitle = paste(
        "Percentage distribution within each",
        "trajectory class"
      ),
      x = NULL,
      y = "Percentage of patients",
      fill = legend_title
    ) +
    characteristics_theme
}



############################################################
# FOREST-PLOT FUNCTION
############################################################

make_forest_plot <- function(
    data,
    class_summary,
    selected_outcome,
    selected_model,
    biomarker
) {
  
  forest_data <- data %>%
    filter(
      outcome == selected_outcome,
      model == selected_model
    )
  
  ##########################################################
  # DEFINE CLASS NAMES AND REFERENCE CLASS
  ##########################################################
  
  if (biomarker == "pp") {
    
    forest_data <- forest_data %>%
      mutate(
        trajectory_label = gsub(
          dashboard_config$biomarker_a$term_prefix,
          "",
          term
        )
      )
    
    reference_class <- dashboard_config$biomarker_a$class_labels[[1]]
    
    class_order <- dashboard_config$biomarker_a$class_labels
    
    class_n_lookup <- setNames(
      class_summary$patients,
      as.character(
        class_summary$trajectory_class
      )
    )
    
  } else {
    
    forest_data <- forest_data %>%
      mutate(
        trajectory_label = gsub(
          dashboard_config$biomarker_b$term_prefix,
          "",
          term
        )
      )
    
    reference_class <- dashboard_config$biomarker_b$class_labels[[1]]
    
    class_order <- dashboard_config$biomarker_b$class_labels
    class_n_lookup <- setNames(
      class_summary$patients,
      as.character(
        class_summary$trajectory_b2m_class
      )
    )
  }
  
  
  ##########################################################
  # ADD THE REFERENCE CLASS
  ##########################################################
  
  reference_row <- tibble(
    outcome = selected_outcome,
    model = selected_model,
    term = NA_character_,
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA_real_,
    trajectory_label = reference_class
  )
  
  
  forest_data <- bind_rows(
    reference_row,
    forest_data
  ) %>%
    mutate(
      is_reference =
        trajectory_label == reference_class,
      
      trajectory_label = factor(
        trajectory_label,
        levels = rev(class_order)
      ),
      
      result_label = ifelse(
        is_reference,
        
        "1.00 (reference)",
        
        paste0(
          sprintf(
            "%.2f",
            estimate
          ),
          " (",
          sprintf(
            "%.2f",
            conf.low
          ),
          "–",
          sprintf(
            "%.2f",
            conf.high
          ),
          "); p=",
          vapply(
            p.value,
            format_p_value,
            character(1)
          )
        )
      )
    )
  
  
  ##########################################################
  # CALCULATE PLOT LIMITS
  ##########################################################
  
  minimum_ci <- min(
    forest_data$conf.low,
    na.rm = TRUE
  )
  
  maximum_ci <- max(
    forest_data$conf.high,
    na.rm = TRUE
  )
  
  lower_limit <- min(
    0.5,
    minimum_ci * 0.8
  )
  
  label_position <- maximum_ci * 1.25
  
  upper_limit <- maximum_ci * 3.2
  
  
  ##########################################################
  # CREATE FOREST PLOT
  ##########################################################
  
  ggplot(
    forest_data,
    aes(
      y = trajectory_label
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#607D8B",
      linewidth = 0.8
    ) +
    
    geom_segment(
      data = forest_data %>%
        filter(
          !is_reference
        ),
      
      aes(
        x = conf.low,
        xend = conf.high,
        yend = trajectory_label
      ),
      
      colour = "#0B2E59",
      linewidth = 1.3
    ) +
    
    geom_point(
      aes(
        x = estimate,
        shape = is_reference
      ),
      
      colour = "#1565C0",
      fill = "#1565C0",
      size = 4
    ) +
    
    geom_text(
      aes(
        x = label_position,
        label = result_label
      ),
      
      hjust = 0,
      size = 4.2,
      colour = "#17212B"
    ) +
    
    scale_shape_manual(
      values = c(
        "FALSE" = 16,
        "TRUE" = 18
      )
    ) +
    scale_y_discrete(
      labels = function(class_name) {
        
        paste0(
          class_name,
          "\n(n = ",
          format(
            class_n_lookup[class_name],
            big.mark = ",",
            scientific = FALSE
          ),
          ")"
        )
      }
    ) +
    
    scale_x_log10(
      limits = c(
        lower_limit,
        upper_limit
      ),
      
      breaks = c(
        0.5,
        0.75,
        1,
        1.5,
        2,
        3,
        4,
        5
      )
    ) +
    
    labs(
      title = paste(
        selected_model,
        selected_outcome
      ),
      
      subtitle = paste(
        "Hazard ratios and 95% confidence intervals;",
        reference_class,
        "is the reference class"
      ),
      
      x = "Hazard ratio, logarithmic scale",
      
      y = NULL
    ) +
    
    guides(
      shape = "none"
    ) +
    
    theme_minimal(
      base_size = 15
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 19,
        colour = "#17212B"
      ),
      
      plot.subtitle = element_text(
        colour = "#607D8B",
        margin = margin(
          b = 20
        )
      ),
      
      axis.text.y = element_text(
        size = 13,
        colour = "#17212B"
      ),
      
      axis.title.x = element_text(
        face = "bold"
      ),
      
      panel.grid.major.y = element_blank(),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        20,
        40,
        20,
        20
      )
    )
}




dashboard_image_folder <- "www"

############################################################
# 7. USER INTERFACE
############################################################

ui <- page_navbar(
  
  title = dashboard_config$app_title,
  theme = dashboard_theme,
  fillable = TRUE,
  
  tags$head(
    tags$style(
      HTML(
        "
      .dashboard-image {
        width: 100%;
        overflow: hidden;
      }

      .dashboard-image img {
        width: 100% !important;
        height: 100% !important;
        object-fit: contain;
      }
      "
      )
    )
  ),
  
  ##########################################################
  # OVERVIEW
  ##########################################################
  
  nav_panel(
    title = "Overview",
    
    div(
      class = "p-4",
      
      h1("Biomarker trajectories and clinical outcomes"),
      
      p(
        class = "lead",
        paste(
          "Explore", dashboard_config$biomarker_a$name, "and",
          dashboard_config$biomarker_b$name,
          "trajectory patterns observed during",
          dashboard_config$trajectory_window, "in",
          paste0(dashboard_config$study_name, ".")
        )
      ),
      
      h3("Full cohort"),
      
      uiOutput(
        outputId = "full_cohort_summary"
      ),
      
      h3(
        class = "mt-4",
        "Cohort characteristics"
      ),
      
      layout_columns(
        
        card(
          full_screen = TRUE,
          card_body(
            div(
              class = "dashboard-image",
              imageOutput(
                outputId = "full_cohort_age_plot",
                width = "100%",
                height = "430px"
              )
            )
          )
        ),
        
        card(
          full_screen = TRUE,
          card_body(
            div(
              class = "dashboard-image",
              imageOutput(
                outputId = "full_cohort_sex_plot",
                width = "100%",
                height = "430px"
              )
            )
          )
        ),
        
        card(
          full_screen = TRUE,
          card_body(
            div(
              class = "dashboard-image",
              imageOutput(
                outputId = "full_cohort_r2iss_plot",
                width = "100%",
                height = "430px"
              )
            )
          )
        ),
        
        card(
          full_screen = TRUE,
          card_body(
            div(
              class = "dashboard-image",
              imageOutput(
                outputId = "full_cohort_pathway_plot",
                width = "100%",
                height = "430px"
              )
            )
          )
        ),
        
        col_widths = c(6, 6, 6, 6)
      ),
      
      layout_columns(
        
        value_box(
          title = "Biomarkers",
          value = "2",
          showcase = bsicons::bs_icon("activity"),
          theme = "primary"
        ),
        
        value_box(
          title = paste(dashboard_config$biomarker_a$short_name, "cohort"),
          value = tags$div(
            style = "font-size: 1.8rem; font-weight: 400;",
            dashboard_config$biomarker_a$cohort_n
          ),
          showcase = bsicons::bs_icon("graph-up"),
          theme = "success"
        ),
        
        value_box(
          title = paste(dashboard_config$biomarker_b$short_name, "cohort"),
          value = tags$div(
            style = "font-size: 1.8rem; font-weight: 400;",
            dashboard_config$biomarker_b$cohort_n
          ),
          showcase = bsicons::bs_icon("diagram-3"),
          theme = "secondary"
        ),
        
        col_widths = c(4, 4, 4)
      ),
      
      layout_columns(
        
        card(
          card_header(paste(dashboard_config$biomarker_a$name, "analysis")),
          card_body(
            h4(paste(
              length(dashboard_config$biomarker_a$class_labels),
              "trajectory classes"
            )),
            p(
              paste(
                "Explore trajectory patterns,",
                "compare patient characteristics across classes,",
                "and examine subsequent survival and adjusted outcomes."
              )
            )
          )
        ),
        
        card(
          card_header(paste(dashboard_config$biomarker_b$name, "analysis")),
          card_body(
            h4(paste(
              length(dashboard_config$biomarker_b$class_labels),
              "trajectory classes"
            )),
            p(
              paste(
                "Explore trajectory patterns, compare patient",
                "characteristics across classes, and examine subsequent",
                "survival and adjusted outcomes."
              )
            )
          )
        ),
        
        col_widths = c(6, 6)
      ),
      card(
        class = "mt-3",
        card_header("Analytical approach"),
        card_body(
          tags$ul(
            tags$li(
              paste(
                "Repeated biomarker measurements during the selected window",
                "were used to identify patterns of biomarker change during treatment."
              )
            ),
            tags$li(
              paste(
                "LCGA identified groups with similar biomarker patterns.",
                "The number of groups was selected using model fit,",
                "group separation, size and clinical interpretation."
              )
            ),
            tags$li(
              paste(
                "Progression and survival after the trajectory window were compared",
                "across the groups using Kaplan–Meier and Cox models."
              )
            )
          ),
          p(
            class = "mb-0 text-muted",
            paste(
              "The groups were identified from this dataset and are not",
              "established clinical treatment categories."
            )
          )
        )
      )
    )
  ),
  
  ##########################################################
  # PARAPROTEIN
  ##########################################################
  
  nav_panel(
    title = dashboard_config$biomarker_a$short_name,
    
    div(
      class = "p-4",
      
      h2(paste(dashboard_config$biomarker_a$name, "trajectory analysis")),
      
      p(
        class = "lead",
        paste(
          length(dashboard_config$biomarker_a$class_labels),
          "trajectory classes identified during",
          paste0(dashboard_config$trajectory_window, ".")
        )
      ),
      
      navset_card_tab(
        
        nav_panel(
          title = "Trajectory patterns",
          
          layout_columns(
            
            card(
              class = "trajectory-card",
              full_screen = TRUE,
              card_header(
                paste(
                  dashboard_config$biomarker_a$name,
                  "trajectories during",
                  dashboard_config$trajectory_window
                )
              ),
              card_body(
                plotlyOutput(
                  outputId = "pp_trajectory_plot",
                  height = "600px"
                )
              ),
              card_footer(
                p(
                  class = "mb-0 text-muted small",
                  paste(
                    "Dots show the observed average biomarker measurements",
                    "at each time point. Lines show the model-estimated",
                    "trajectory for each patient group."
                  )
                )
              )
            ),
            
            uiOutput(
              outputId = "pp_selected_class_summary"
            ),
            
            col_widths = c(8, 4)
          )
        ),
        
        nav_panel(
          title = "Class characteristics",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 300,
              h4("Select a characteristic"),
              selectInput(
                inputId = "pp_characteristic",
                label = NULL,
                choices = setNames(
                  c(
                    "secretory", "light_chain", "r2iss",
                    "pathway", "sex", "age"
                  ),
                  c(
                    dashboard_config$secretory_type_label,
                    dashboard_config$light_chain_label,
                    dashboard_config$stage_label,
                    dashboard_config$pathway_label,
                    "Sex", "Age at baseline"
                  )
                ),
                selected = "secretory"
              ),
              hr(),
              p(
                class = "text-muted",
                paste(
                  "Select a characteristic to compare its distribution",
                  "across the trajectory classes."
                )
              )
            ),
            
            card(
              full_screen = TRUE,
              card_body(
                plotOutput(
                  outputId = "pp_characteristics_plot",
                  height = "600px"
                )
              )
            )
          )
        ),
        
        nav_panel(
          title = "Survival",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 260,
              radioButtons(
                inputId = "pp_survival_outcome",
                label = "Select outcome",
                choices = c(
                  "Progression-free survival" = "pfs",
                  "Overall survival" = "os"
                ),
                selected = "pfs"
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header(
                textOutput(
                  outputId = "pp_survival_title"
                )
              ),
              card_body(
                div(
                  class = "dashboard-image",
                  
                  imageOutput(
                    outputId = "pp_survival_plot",
                    width = "100%",
                    height = "650px"
                  )
                )
              )
            )
          )
        ),
        
        nav_panel(
          title = "Adjusted analysis",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 300,
              radioButtons(
                inputId = "pp_forest_outcome",
                label = "Select outcome",
                choices = c(
                  "Progression-free survival" =
                    "Progression-free survival",
                  "Overall survival" =
                    "Overall survival"
                ),
                selected = "Progression-free survival"
              ),
              hr(),
              radioButtons(
                inputId = "pp_forest_model",
                label = "Select model",
                choices = c(
                  "Unadjusted" = "Unadjusted",
                  "Adjusted" = "Adjusted"
                ),
                selected = "Adjusted"
              ),
              hr(),
              p(
                class = "text-muted",
                paste(
                  "Adjusted models include",
                  paste0(dashboard_config$adjusted_covariates, ".")
                )
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header(paste(
                dashboard_config$biomarker_a$name,
                "trajectory hazard ratios"
              )),
              card_body(
                plotOutput(
                  outputId = "pp_forest_plot",
                  height = "600px"
                )
              )
            )
          )
        )
      )
    )
  ),
  
  ##########################################################
  # B2M
  ##########################################################
  
  nav_panel(
    title = dashboard_config$biomarker_b$short_name,
    
    div(
      class = "p-4",
      
      h2(paste(dashboard_config$biomarker_b$name, "trajectory analysis")),
      
      p(
        class = "lead",
        paste(
          length(dashboard_config$biomarker_b$class_labels),
          "trajectory classes identified during",
          paste0(dashboard_config$trajectory_window, ".")
        )
      ),
      
      navset_card_tab(
        
        nav_panel(
          title = "Trajectory patterns",
          
          layout_columns(
            
            card(
              class = "trajectory-card",
              full_screen = TRUE,
              card_header(
                paste(
                  dashboard_config$biomarker_b$name,
                  "trajectories during",
                  dashboard_config$trajectory_window
                )
              ),
              card_body(
                plotlyOutput(
                  outputId = "b2m_trajectory_plot",
                  height = "600px"
                )
              ),
              card_footer(
                p(
                  class = "mb-0 text-muted small",
                  paste(
                    "Dots show the observed average biomarker",
                    "measurements at each time point. Lines show the",
                    "model-estimated trajectory for each patient group."
                  )
                )
              )
            ),
            
            uiOutput(
              outputId = "b2m_selected_class_summary"
            ),
            
            col_widths = c(8, 4)
          )
        ),
        
        nav_panel(
          title = "Class characteristics",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 300,
              h4("Select a characteristic"),
              selectInput(
                inputId = "b2m_characteristic",
                label = NULL,
                choices = setNames(
                  c("r2iss", "pathway", "sex", "age"),
                  c(
                    dashboard_config$stage_label,
                    dashboard_config$pathway_label,
                    "Sex", "Age at baseline"
                  )
                ),
                selected = "r2iss"
              ),
              hr(),
              p(
                class = "text-muted",
                paste(
                  "Select a characteristic to compare its distribution",
                  "across the trajectory classes."
                )
              )
            ),
            
            card(
              full_screen = TRUE,
              card_body(
                plotOutput(
                  outputId = "b2m_characteristics_plot",
                  height = "600px"
                )
              ),
              card_footer(
                p(
                  class = "mb-0 text-muted",
                  paste(
                    "Percentages are calculated within each",
                    "trajectory class."
                  )
                )
              )
            )
          )
        ),
        
        nav_panel(
          title = "Survival",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 260,
              radioButtons(
                inputId = "b2m_survival_outcome",
                label = "Select outcome",
                choices = c(
                  "Progression-free survival" = "pfs",
                  "Overall survival" = "os"
                ),
                selected = "pfs"
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header(
                textOutput(
                  outputId = "b2m_survival_title"
                )
              ),
              card_body(
                div(
                  class = "dashboard-image",
                  
                  imageOutput(
                    outputId = "b2m_survival_plot",
                    width = "100%",
                    height = "650px"
                  )
                )
              )
            )
          )
        ),
        
        nav_panel(
          title = "Adjusted analysis",
          
          layout_sidebar(
            sidebar = sidebar(
              width = 300,
              radioButtons(
                inputId = "b2m_forest_outcome",
                label = "Select outcome",
                choices = c(
                  "Progression-free survival" =
                    "Progression-free survival",
                  "Overall survival" =
                    "Overall survival"
                ),
                selected = "Progression-free survival"
              ),
              hr(),
              radioButtons(
                inputId = "b2m_forest_model",
                label = "Select model",
                choices = c(
                  "Unadjusted" = "Unadjusted",
                  "Adjusted" = "Adjusted"
                ),
                selected = "Adjusted"
              ),
              hr(),
              p(
                class = "text-muted",
                paste(
                  "Adjusted models include",
                  paste0(dashboard_config$adjusted_covariates, ".")
                )
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header(
                paste(
                  dashboard_config$biomarker_b$name,
                  "trajectory hazard ratios"
                )
              ),
              card_body(
                plotOutput(
                  outputId = "b2m_forest_plot",
                  height = "600px"
                )
              )
            )
          )
        )
      )
    )
  )
)


############################################################
# 8. SERVER
############################################################

server <- function(input, output, session) {
  
  ##########################################################
  # HELPER: ADD CLASS NAME TO EVERY PLOTLY TRACE
  ##########################################################
  
  make_clickable_trajectory <- function(plot_object, plot_source) {
    
    interactive_plot <- ggplotly(
      plot_object,
      tooltip = c("colour", "x", "y"),
      source = plot_source
    )
    
    # to attach the trace/class name to every plotted point.
    # event_data() can then identify the class that was clicked.
    for (trace_number in seq_along(interactive_plot$x$data)) {
      
      trace_name <- interactive_plot$x$data[[trace_number]]$name
      number_of_points <- length(
        interactive_plot$x$data[[trace_number]]$x
      )
      
      interactive_plot$x$data[[trace_number]]$customdata <- rep(
        trace_name,
        number_of_points
      )
    }
    
    interactive_plot <- event_register(
      interactive_plot,
      "plotly_click"
    )
    interactive_plot %>%
      config(
        displayModeBar = TRUE,
        responsive = TRUE
      ) %>%
      layout(
        autosize = TRUE
      )
  }
  
  ##########################################################
  # HELPER: MATCH A CLICKED TRACE TO AN ALLOWED CLASS
  ##########################################################
  
  identify_clicked_class <- function(clicked_name, allowed_classes) {
    
    if (is.null(clicked_name) || is.na(clicked_name)) {
      return(NULL)
    }
    
    matched_class <- allowed_classes[
      vapply(
        allowed_classes,
        function(class_name) {
          grepl(
            class_name,
            clicked_name,
            fixed = TRUE
          )
        },
        logical(1)
      )
    ]
    
    if (length(matched_class) == 0) {
      return(NULL)
    }
    
    matched_class[[1]]
  }
  
  ##########################################################
  # FULL-COHORT SUMMARY
  ##########################################################
  
  output$full_cohort_summary <- renderUI({
    
    cohort_summary <- dashboard_outputs$full_cohort$summary
    
    req(
      nrow(cohort_summary) == 1
    )
    
    r2iss_percentage <-
      cohort_summary$r2iss_available /
      cohort_summary$patients * 100
    
    layout_columns(
      
      value_box(
        title = "Full cohort",
        value = paste(
          format(
            cohort_summary$patients,
            big.mark = ","
          ),
          "patients"
        ),
        theme = "primary"
      ),
      
      value_box(
        title = "Median age",
        value = paste0(
          round(cohort_summary$median_age),
          " years"
        ),
        p(
          paste0(
            "Q1–Q3: ",
            round(cohort_summary$q1_age),
            "–",
            round(cohort_summary$q3_age)
          )
        ),
        theme = "secondary"
      ),
      
      value_box(
        title = paste(dashboard_config$stage_label, "recorded"),
        value = paste(
          format(
            cohort_summary$r2iss_available,
            big.mark = ","
          ),
          "patients"
        ),
        p(
          paste0(
            round(r2iss_percentage, 1),
            "% of cohort"
          )
        ),
        theme = "success"
      ),
      
      col_widths = c(4, 4, 4)
    )
  })
  
  output$full_cohort_age_plot <- renderImage({
    
    list(
      src = file.path(
        dashboard_image_folder,
        dashboard_config$files$full_cohort_age
      ),
      contentType = "image/png",
      alt = "Age distribution in the full cohort"
    )
    
  }, deleteFile = FALSE)
  
  
  output$full_cohort_sex_plot <- renderImage({
    
    list(
      src = file.path(
        dashboard_image_folder,
        dashboard_config$files$full_cohort_sex
      ),
      contentType = "image/png",
      alt = "Sex distribution in the full cohort"
    )
    
  }, deleteFile = FALSE)
  
  
  output$full_cohort_r2iss_plot <- renderImage({
    
    list(
      src = file.path(
        dashboard_image_folder,
        dashboard_config$files$full_cohort_stage
      ),
      contentType = "image/png",
      alt = paste(dashboard_config$stage_label, "distribution in the full cohort")
    )
    
  }, deleteFile = FALSE)
  
  
  output$full_cohort_pathway_plot <- renderImage({
    
    list(
      src = file.path(
        dashboard_image_folder,
        dashboard_config$files$full_cohort_pathway
      ),
      contentType = "image/png",
      alt = "Treatment pathway distribution in the full cohort"
    )
    
  }, deleteFile = FALSE)
  
  ##########################################################
  # PP TRAJECTORY AND CLICK SELECTION
  ##########################################################
  
  pp_class_levels <- dashboard_config$biomarker_a$class_labels
  
  selected_pp_class <- reactiveVal(
    dashboard_config$biomarker_a$class_labels[[1]]
  )
  
  output$pp_trajectory_plot <- renderPlotly({
    
    make_clickable_trajectory(
      plot_object = trajectory_plots$pp,
      plot_source = "pp_plot"
    )
  })
  
  observeEvent(
    event_data(
      "plotly_click",
      source = "pp_plot"
    ),
    {
      click_information <- event_data(
        "plotly_click",
        source = "pp_plot"
      )
      
      req(
        nrow(click_information) >= 1
      )
      
      clicked_class <- identify_clicked_class(
        clicked_name = as.character(
          click_information$customdata[[1]]
        ),
        allowed_classes = pp_class_levels
      )
      
      if (!is.null(clicked_class)) {
        selected_pp_class(clicked_class)
      }
    },
    ignoreInit = TRUE
  )
  
  output$pp_selected_class_summary <- renderUI({
    
    selected_summary <- dashboard_outputs$pp$class_summary %>%
      filter(
        as.character(trajectory_class) ==
          selected_pp_class()
      )
    
    req(
      nrow(selected_summary) == 1
    )
    
    card(
      class = "h-100",
      card_header(
        paste0(
          selected_pp_class(),
          " characteristics"
        )
      ),
      card_body(
        p(
          class = "text-muted",
          "Click a trajectory line or point to select another class."
        ),
        value_box(
          title = "Patients",
          value = format(
            selected_summary$patients,
            big.mark = ","
          ),
          theme = "primary"
        ),
        value_box(
          title = "Median age",
          value = paste0(
            round(selected_summary$median_age),
            " years"
          ),
          p(
            paste0(
              "Q1–Q3: ",
              round(selected_summary$q1_age),
              "–",
              round(selected_summary$q3_age),
              " years"
            )
          ),
          theme = "secondary"
        ),
        value_box(
          title = "Male",
          value = paste0(
            round(selected_summary$male_percent, 1),
            "%"
          ),
          theme = "success"
        ),
        value_box(
          title = "Transplant eligible",
          value = paste0(
            round(
              selected_summary$transplant_eligible_percent,
              1
            ),
            "%"
          ),
          theme = "primary"
        )
      )
    )
  })
  
  ##########################################################
  # PP CHARACTERISTICS
  ##########################################################
  
  output$pp_characteristics_plot <- renderPlot({
    
    req(
      input$pp_characteristic
    )
    
    make_characteristics_plot(
      summary_data = pp_characteristics_summary,
      class_column = "trajectory_class",
      class_levels = pp_class_levels,
      characteristic = input$pp_characteristic
    )
  },
  height = 600,
  res = 110)
  
  ##########################################################
  # PP SURVIVAL
  ##########################################################
  
  output$pp_survival_title <- renderText({
    
    if (input$pp_survival_outcome == "pfs") {
      paste(dashboard_config$landmark_label, "progression-free survival")
    } else {
      paste(dashboard_config$landmark_label, "overall survival")
    }
  })
  
  output$pp_survival_plot <- renderImage({
    
    filename <- if (
      input$pp_survival_outcome == "pfs"
    ) {
      dashboard_config$files$biomarker_a_pfs
    } else {
      dashboard_config$files$biomarker_a_os
    }
    
    list(
      src = file.path(
        dashboard_image_folder,
        filename
      ),
      contentType = "image/png",
      alt = paste(dashboard_config$biomarker_a$name, "landmark survival curve")
    )
    
  }, deleteFile = FALSE)
  
  ##########################################################
  # PP COX FOREST PLOT
  ##########################################################
  
  output$pp_forest_plot <- renderPlot({
    
    req(
      input$pp_forest_outcome,
      input$pp_forest_model
    )
    
    make_forest_plot(
      data = dashboard_outputs$pp$cox_results,
      class_summary = dashboard_outputs$pp$class_summary,
      selected_outcome = input$pp_forest_outcome,
      selected_model = input$pp_forest_model,
      biomarker = "pp"
    )
  },
  height = 600,
  res = 110)
  
  ##########################################################
  # B2M TRAJECTORY AND CLICK SELECTION
  ##########################################################
  
  b2m_class_levels <- dashboard_config$biomarker_b$class_labels
  
  selected_b2m_class <- reactiveVal(
    dashboard_config$biomarker_b$class_labels[[1]]
  )
  
  output$b2m_trajectory_plot <- renderPlotly({
    
    make_clickable_trajectory(
      plot_object = trajectory_plots$b2m,
      plot_source = "b2m_plot"
    )
  })
  
  observeEvent(
    event_data(
      "plotly_click",
      source = "b2m_plot"
    ),
    {
      click_information <- event_data(
        "plotly_click",
        source = "b2m_plot"
      )
      
      req(
        nrow(click_information) >= 1
      )
      
      clicked_class <- identify_clicked_class(
        clicked_name = as.character(
          click_information$customdata[[1]]
        ),
        allowed_classes = b2m_class_levels
      )
      
      if (!is.null(clicked_class)) {
        selected_b2m_class(clicked_class)
      }
    },
    ignoreInit = TRUE
  )
  
  output$b2m_selected_class_summary <- renderUI({
    
    selected_summary <- dashboard_outputs$b2m$class_summary %>%
      filter(
        as.character(trajectory_b2m_class) ==
          selected_b2m_class()
      )
    
    req(
      nrow(selected_summary) == 1
    )
    
    card(
      class = "h-100",
      card_header(
        paste0(
          selected_b2m_class(),
          " characteristics"
        )
      ),
      card_body(
        p(
          class = "text-muted",
          "Click a trajectory line or point to select another class."
        ),
        value_box(
          title = "Patients",
          value = format(
            selected_summary$patients,
            big.mark = ","
          ),
          theme = "primary"
        ),
        value_box(
          title = "Median age",
          value = paste0(
            round(selected_summary$median_age),
            " years"
          ),
          p(
            paste0(
              "Q1–Q3: ",
              round(selected_summary$q1_age),
              "–",
              round(selected_summary$q3_age),
              " years"
            )
          ),
          theme = "secondary"
        ),
        value_box(
          title = "Male",
          value = paste0(
            round(selected_summary$male_percent, 1),
            "%"
          ),
          theme = "success"
        ),
        value_box(
          title = "Transplant eligible",
          value = paste0(
            round(
              selected_summary$transplant_eligible_percent,
              1
            ),
            "%"
          ),
          theme = "primary"
        )
      )
    )
  })
  
  ##########################################################
  # B2M CHARACTERISTICS
  ##########################################################
  
  output$b2m_characteristics_plot <- renderPlot({
    
    req(
      input$b2m_characteristic
    )
    
    make_characteristics_plot(
      summary_data = b2m_characteristics_summary,
      class_column = "trajectory_b2m_class",
      class_levels = b2m_class_levels,
      characteristic = input$b2m_characteristic
    )
  },
  height = 600,
  res = 110)
  
  ##########################################################
  # B2M SURVIVAL
  ##########################################################
  
  output$b2m_survival_title <- renderText({
    
    if (input$b2m_survival_outcome == "pfs") {
      paste(dashboard_config$landmark_label, "progression-free survival")
    } else {
      paste(dashboard_config$landmark_label, "overall survival")
    }
  })
  
  output$b2m_survival_plot <- renderImage({
    
    filename <- if (
      input$b2m_survival_outcome == "pfs"
    ) {
      dashboard_config$files$biomarker_b_pfs
    } else {
      dashboard_config$files$biomarker_b_os
    }
    
    list(
      src = file.path(
        dashboard_image_folder,
        filename
      ),
      contentType = "image/png",
      alt = paste(dashboard_config$biomarker_b$name, "landmark survival curve")
    )
    
  }, deleteFile = FALSE)
  
  ##########################################################
  # B2M COX FOREST PLOT
  ##########################################################
  
  output$b2m_forest_plot <- renderPlot({
    
    req(
      input$b2m_forest_outcome,
      input$b2m_forest_model
    )
    
    make_forest_plot(
      data = dashboard_outputs$b2m$cox_results,
      class_summary = dashboard_outputs$b2m$class_summary,
      selected_outcome = input$b2m_forest_outcome,
      selected_model = input$b2m_forest_model,
      biomarker = "b2m"
    )
  },
  height = 600,
  res = 110)
}

############################################################
# 9. RUN APPLICATION
############################################################

shinyApp(
  ui = ui,
  server = server
)
