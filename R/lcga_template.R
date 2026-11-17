###############################################################################
# REUSABLE LATENT CLASS GROWTH ANALYSIS TEMPLATE
###############################################################################
# SUMMARY
# This script identifies groups whose repeated biomarker measurements follow
# similar patterns over time. It:
#   1. fits quadratic LCGA models containing one to four classes;
#   2. compares the candidate models using BIC, entropy and class sizes;
#   3. assigns each patient to the class with the highest posterior probability;
#   4. creates a publication-ready trajectory plot.
#
# SETUP
# Edit only Section 2. The script maps the user's column names to stable
# internal names before fitting the models.
#
# DATA GOVERNANCE
# The source data are controlled-access clinical-trial data and are not included
# in the GitHub repository. This script does not write patient-level data or
# model objects to disk.
#
# MAIN OUTPUTS
#   - model_selection: numerical comparison of the candidate LCGA models;
#   - model_selection_table: formatted model-selection table;
#   - trajectory_classification: model-derived class assignment in memory;
#   - biomarker_trajectory_plot: predicted trajectories and weighted observed
#     means.
###############################################################################


# 1. PACKAGES -----------------------------------------------------------------

library(lcmm)
library(tidyverse)
library(gt)


# 2. TEMPLATE SETTINGS — EDIT THIS SECTION ONLY -------------------------------

lcga_config <- list(
  input_object = "YOUR_PREPARED_DATA_OBJECT",
  patient_id_column = "YOUR_PATIENT_ID_COLUMN",
  time_column = "YOUR_TIME_COLUMN",
  biomarker_column = "YOUR_BIOMARKER_COLUMN",
  biomarker_name = "YOUR BIOMARKER FULL NAME",
  biomarker_short_name = "YOUR BIOMARKER",
  biomarker_unit = "YOUR UNIT",
  time_window_label = "YOUR ANALYSIS WINDOW",
  time_min = 0,
  time_max = 3,
  selected_classes = 3,
  class_labels = c(
    "CLASS 1 LABEL",
    "CLASS 2 LABEL",
    "CLASS 3 LABEL"
  ),
  class_colours = c(
    "#1565C0",
    "#00897B",
    "#0B2E59"
  ),
  random_seed = 123
)

# The named object must already exist in the R environment. Only prepared,
# approved data should be used; no source data are included in this repository.
if (!exists(lcga_config$input_object, envir = .GlobalEnv)) {
  stop(
    "Create or load the prepared data object named: ",
    lcga_config$input_object
  )
}

input_data <- get(
  lcga_config$input_object,
  envir = .GlobalEnv
)

required_columns <- c(
  lcga_config$patient_id_column,
  lcga_config$time_column,
  lcga_config$biomarker_column
)

missing_columns <- setdiff(required_columns, names(input_data))

if (length(missing_columns) > 0) {
  stop(
    "The input data are missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (
  length(lcga_config$class_labels) != lcga_config$selected_classes ||
  length(lcga_config$class_colours) != lcga_config$selected_classes
) {
  stop(
    "Provide one class label and one colour for every selected class."
  )
}

if (!lcga_config$selected_classes %in% 1:4) {
  stop("selected_classes must be one of 1, 2, 3 or 4.")
}

# Map user-supplied columns to the stable internal names used below.
lcga_b2m <- input_data %>%
  transmute(
    pat_no = .data[[lcga_config$patient_id_column]],
    sample_datem = as.numeric(.data[[lcga_config$time_column]]),
    b2m = as.numeric(.data[[lcga_config$biomarker_column]])
  ) %>%
  filter(
    !is.na(pat_no),
    !is.na(sample_datem),
    !is.na(b2m),
    sample_datem >= lcga_config$time_min,
    sample_datem <= lcga_config$time_max
  ) %>%
  as.data.frame()


# 3. FIT CANDIDATE LCGA MODELS -------------------------------------------------
# All models use a quadratic function of time and assume no within-class random
# effects, consistent with latent class growth analysis.

lcga_b2m_1_class <- hlme(
  fixed = b2m ~ sample_datem + I(sample_datem^2),
  random = ~ -1,
  subject = "pat_no",
  ng = 1,
  data = lcga_b2m,
  na.action = 1
)

set.seed(lcga_config$random_seed)
lcga_b2m_2_class <- hlme(
  fixed = b2m ~ sample_datem + I(sample_datem^2),
  mixture = ~ sample_datem + I(sample_datem^2),
  random = ~ -1,
  subject = "pat_no",
  ng = 2,
  data = lcga_b2m,
  B = lcga_b2m_1_class,
  na.action = 1
)

set.seed(lcga_config$random_seed)
lcga_b2m_3_class <- hlme(
  fixed = b2m ~ sample_datem + I(sample_datem^2),
  mixture = ~ sample_datem + I(sample_datem^2),
  random = ~ -1,
  subject = "pat_no",
  ng = 3,
  data = lcga_b2m,
  B = lcga_b2m_1_class,
  na.action = 1
)

set.seed(lcga_config$random_seed)
lcga_b2m_4_class <- hlme(
  fixed = b2m ~ sample_datem + I(sample_datem^2),
  mixture = ~ sample_datem + I(sample_datem^2),
  random = ~ -1,
  subject = "pat_no",
  ng = 4,
  data = lcga_b2m,
  B = lcga_b2m_1_class,
  na.action = 1
)

candidate_models <- list(
  lcga_b2m_1_class,
  lcga_b2m_2_class,
  lcga_b2m_3_class,
  lcga_b2m_4_class
)

# Inspect convergence codes. A value of 1 indicates normal convergence in
# `lcmm`; all candidate models should be checked before interpretation.
convergence_codes <- tibble(
  classes = 1:4,
  convergence_code = map_int(candidate_models, ~ .x$conv)
)

print(convergence_codes)


# 4. COMPARE THE CANDIDATE MODELS ---------------------------------------------

calculate_entropy <- function(model) {
  number_of_classes <- model$ng
  
  if (number_of_classes == 1) {
    return(NA_real_)
  }
  
  probabilities <- as.matrix(
    model$pprob[paste0("prob", seq_len(number_of_classes))]
  )
  
  positive_probabilities <- probabilities[probabilities > 0]
  
  1 + sum(
    positive_probabilities * log(positive_probabilities)
  ) / (
    nrow(model$pprob) * log(number_of_classes)
  )
}

format_class_sizes <- function(model) {
  number_of_classes <- model$ng
  
  class_counts <- table(
    factor(
      model$pprob$class,
      levels = seq_len(number_of_classes)
    )
  )
  
  class_counts <- as.numeric(class_counts)
  class_percentages <- class_counts / sum(class_counts) * 100
  
  paste0(
    "Class ", seq_len(number_of_classes), ": ",
    class_counts, " (", round(class_percentages, 1), "%)",
    collapse = "<br>"
  )
}

model_selection <- tibble(
  classes = 1:4,
  BIC = map_dbl(candidate_models, ~ .x$BIC),
  entropy = map_dbl(candidate_models, calculate_entropy),
  class_size = map_chr(candidate_models, format_class_sizes)
)

model_selection_table <- model_selection %>%
  gt() %>%
  cols_label(
    classes = "Number of classes",
    BIC = "BIC",
    entropy = "Entropy",
    class_size = "Class size, n (%)"
  ) %>%
  fmt_number(
    columns = BIC,
    decimals = 2
  ) %>%
  fmt_number(
    columns = entropy,
    decimals = 3
  ) %>%
  sub_missing(
    columns = entropy,
    missing_text = "–"
  ) %>%
  fmt(
    columns = class_size,
    fns = function(value) lapply(value, gt::html)
  ) %>%
  tab_header(
    title = paste(
      lcga_config$biomarker_short_name,
      "LCGA model selection"
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = classes == lcga_config$selected_classes
    )
  ) %>%
  tab_source_note(
    source_note = md(
      paste(
        paste0(
          "The ", lcga_config$selected_classes,
          "-class solution was selected using statistical fit,"
        ),
        "classification quality, class size and clinical interpretability."
      )
    )
  )

model_selection_table


# 5. ASSIGN THE FINAL TRAJECTORY CLASSES --------------------------------------
# The selected model is taken from the candidate-model list. Each patient is
# assigned to the class with the highest posterior probability.

final_model <- candidate_models[[lcga_config$selected_classes]]

class_levels <- lcga_config$class_labels

class_colours <- setNames(
  lcga_config$class_colours,
  class_levels
)

probability_columns <- paste0(
  "prob",
  seq_len(lcga_config$selected_classes)
)

trajectory_classification <- final_model$pprob %>%
  mutate(
    posterior_probability = apply(
      select(., all_of(probability_columns)),
      1,
      max,
      na.rm = TRUE
    )
  ) %>%
  transmute(
    pat_no,
    trajectory_b2m_class = factor(
      class,
      levels = seq_len(lcga_config$selected_classes),
      labels = class_levels
    ),
    posterior_probability
  )

class_distribution <- trajectory_classification %>%
  count(trajectory_b2m_class, name = "patients") %>%
  mutate(
    percentage = patients / sum(patients) * 100
  )

print(class_distribution)


# 6. PREPARE THE TRAJECTORY PLOT ----------------------------------------------

# Obtain the model-predicted trajectory for each class at closely spaced time
# points across the configured analysis window.
prediction_times <- data.frame(
  sample_datem = seq(
    lcga_config$time_min,
    lcga_config$time_max,
    by = 0.01
  )
)

lcga_predictions <- predictY(
  final_model,
  newdata = prediction_times,
  var.time = "sample_datem",
  draws = FALSE
)

prediction_data <- bind_cols(
  prediction_times,
  as.data.frame(lcga_predictions$pred)
) %>%
  pivot_longer(
    cols = starts_with("Ypred_class"),
    names_to = "class_number",
    values_to = "predicted_b2m"
  ) %>%
  mutate(
    class_number = as.integer(
      sub("Ypred_class", "", class_number)
    ),
    trajectory_b2m_class = factor(
      class_levels[class_number],
      levels = class_levels
    )
  )

# Reshape the posterior probabilities so each measurement contributes to the
# observed class mean in proportion to its probability of class membership.
posterior_data <- final_model$pprob %>%
  select(
    pat_no,
    all_of(probability_columns)
  ) %>%
  pivot_longer(
    cols = starts_with("prob"),
    names_to = "class_number",
    values_to = "posterior_probability"
  ) %>%
  mutate(
    class_number = as.integer(
      sub("prob", "", class_number)
    )
  )

# Calculate posterior-probability-weighted observed means in quarter-month
# intervals. These points show the observed data underlying the fitted curves.
time_breaks <- seq(
  lcga_config$time_min,
  lcga_config$time_max,
  by = 0.25
)

observed_data <- lcga_b2m %>%
  left_join(
    posterior_data,
    by = "pat_no"
  ) %>%
  mutate(
    time_interval = cut(
      sample_datem,
      breaks = time_breaks,
      include.lowest = TRUE,
      labels = FALSE
    ),
    time = (
      time_breaks[time_interval] +
        time_breaks[time_interval + 1]
    ) / 2
  ) %>%
  filter(!is.na(time_interval)) %>%
  group_by(class_number, time) %>%
  summarise(
    observed_b2m = weighted.mean(
      b2m,
      w = posterior_probability,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    trajectory_b2m_class = factor(
      class_levels[class_number],
      levels = class_levels
    )
  )

biomarker_trajectory_plot <- ggplot() +
  geom_line(
    data = prediction_data,
    aes(
      x = sample_datem,
      y = predicted_b2m,
      colour = trajectory_b2m_class,
      group = trajectory_b2m_class
    ),
    linewidth = 1.3
  ) +
  geom_point(
    data = observed_data,
    aes(
      x = time,
      y = observed_b2m,
      colour = trajectory_b2m_class
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_colour_manual(values = class_colours) +
  scale_x_continuous(
    breaks = seq(
      lcga_config$time_min,
      lcga_config$time_max,
      by = 0.5
    ),
    limits = c(
      lcga_config$time_min,
      lcga_config$time_max
    ),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = paste(
      lcga_config$biomarker_name,
      "trajectories during",
      lcga_config$time_window_label
    ),
    subtitle = paste(
      "Lines show model-predicted trajectories;",
      "dots show weighted observed means"
    ),
    x = "Time from baseline (months)",
    y = paste0(
      lcga_config$biomarker_name,
      " (", lcga_config$biomarker_unit, ")"
    ),
    colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

biomarker_trajectory_plot

###############################################################################
# END OF SCRIPT
###############################################################################
