###############################################################################
# SYNTHETIC DEMO DATA FOR THE DASHBOARD
###############################################################################
# Creates FAKE aggregate inputs so the dashboard can be run by anyone.
# Every number here is simulated. Nothing comes from any real study.
#
# Run from the repository root:
#   Rscript demo/make_demo_dashboard_data.R
#   shiny::runApp("dashboard")
#
# Output (ignored by git, so real outputs placed here are never committed):
#   dashboard/dashboard_data/*.rds
#   dashboard/www/*.png
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)

set.seed(2026)

data_dir  <- file.path("dashboard", "dashboard_data")
image_dir <- file.path("dashboard", "www")
dir.create(data_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)

# Labels match the placeholder values in dashboard/app.R, so the app runs
# without edits. If you rename classes in app.R, rename them here too.
a_labels <- paste("BIOMARKER A CLASS", 1:4)
b_labels <- paste("BIOMARKER B CLASS", 1:3)
a_prefix <- "trajectory_class"
b_prefix <- "trajectory_b2m_class"

blues <- c("#0B2E59", "#1565C0", "#56A8E0", "#9ED2F2", "#B0B7C3")


# 1. SIMULATE A FAKE PATIENT-LEVEL COHORT -------------------------------------
# Used only to derive consistent fake summaries. It is never saved.

n <- 2000

cohort <- tibble(
  id         = seq_len(n),
  age        = round(rnorm(n, 66, 9)),
  sex        = sample(c("Female", "Male"), n, TRUE, c(0.4, 0.6)),
  r2iss      = sample(c("I", "II", "III", "IV", "Missing"), n, TRUE,
                      c(0.17, 0.32, 0.38, 0.05, 0.08)),
  pathway    = if_else(age < 68, "Transplant Eligible",
                       "Transplant Not Eligible"),
  secretory  = sample(c("Intact immunoglobulin", "Light-chain-only",
                        "Non-secretor"), n, TRUE, c(0.8, 0.17, 0.03)),
  light_chain = sample(c("Kappa", "Lambda", "Not recorded"), n, TRUE,
                       c(0.6, 0.35, 0.05)),
  class_a    = factor(sample(a_labels, n, TRUE, c(0.35, 0.15, 0.45, 0.05)),
                      levels = a_labels),
  class_b    = factor(sample(b_labels, n, TRUE, c(0.5, 0.35, 0.15)),
                      levels = b_labels)
)

# Fake survival times, with mild differences between classes.
hazard <- 0.02 * exp(0.02 * (cohort$age - 66) +
                     0.25 * (as.integer(cohort$class_a) - 1) / 3 +
                     0.20 * (as.integer(cohort$class_b) - 1) / 2)
cohort <- cohort %>%
  mutate(
    pfs_time  = pmin(rexp(n, hazard * 1.8), 120),
    pfs_event = as.integer(pfs_time < 120 & runif(n) < 0.85),
    os_time   = pmin(rexp(n, hazard), 150),
    os_event  = as.integer(os_time < 150 & runif(n) < 0.8)
  )


# 2. TRAJECTORY PLOTS ---------------------------------------------------------

make_trajectory_plot <- function(labels, start, slope, title_unit) {
  times <- seq(0, 3, by = 0.05)
  lines <- expand_grid(class = factor(labels, levels = labels), time = times) %>%
    mutate(k = as.integer(class),
           value = pmax(start[k] * exp(-slope[k] * time), 0.2))
  points <- expand_grid(class = factor(labels, levels = labels),
                        time = seq(0.125, 3, by = 0.25)) %>%
    mutate(k = as.integer(class),
           value = pmax(start[k] * exp(-slope[k] * time), 0.2) *
             exp(rnorm(n(), 0, 0.05)))

  ggplot() +
    geom_line(data = lines,
              aes(time, value, colour = class, group = class),
              linewidth = 1.3) +
    geom_point(data = points,
               aes(time, value, colour = class),
               size = 2, alpha = 0.7) +
    scale_colour_manual(values = setNames(blues[seq_along(labels)], labels)) +
    labs(x = "Time from baseline (months)", y = title_unit, colour = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
}

trajectory_plots <- list(
  pp  = make_trajectory_plot(a_labels, c(8, 25, 40, 30), c(1.5, 0.6, 1.2, 0.1),
                             "Biomarker A (units)"),
  b2m = make_trajectory_plot(b_labels, c(3, 5, 9), c(0.3, 0.4, 0.2),
                             "Biomarker B (units)")
)
saveRDS(trajectory_plots, file.path(data_dir, "trajectory_plots.rds"))


# 3. CLASS CHARACTERISTICS ----------------------------------------------------

summarise_characteristics <- function(class_var, class_name, variables) {
  out <- list()

  out$age <- cohort %>%
    group_by(.data[[class_var]]) %>%
    summarise(minimum = quantile(age, 0.05), lower = quantile(age, 0.25),
              middle = median(age), upper = quantile(age, 0.75),
              maximum = quantile(age, 0.95), .groups = "drop") %>%
    rename(!!class_name := all_of(class_var))

  for (v in variables) {
    out[[v]] <- cohort %>%
      count(.data[[class_var]], category_display = .data[[v]]) %>%
      group_by(.data[[class_var]]) %>%
      mutate(percentage = n / sum(n) * 100) %>%
      ungroup() %>%
      select(-n) %>%
      rename(!!class_name := all_of(class_var))
  }
  out
}

all_vars <- c("sex", "r2iss", "pathway", "secretory", "light_chain")

saveRDS(summarise_characteristics("class_a", "trajectory_class", all_vars),
        file.path(data_dir, "biomarker_a_characteristics.rds"))
saveRDS(summarise_characteristics("class_b", "trajectory_b2m_class",
                                  c("sex", "r2iss", "pathway")),
        file.path(data_dir, "biomarker_b_characteristics.rds"))


# 4. CLASS SUMMARIES AND COX RESULTS ------------------------------------------

class_summary <- function(class_var, class_name) {
  cohort %>%
    group_by(.data[[class_var]]) %>%
    summarise(patients = n(), median_age = median(age),
              q1_age = quantile(age, 0.25), q3_age = quantile(age, 0.75),
              male_percent = mean(sex == "Male") * 100,
              transplant_eligible_percent =
                mean(pathway == "Transplant Eligible") * 100,
              .groups = "drop") %>%
    rename(!!class_name := all_of(class_var))
}

cox_results <- function(class_var, prefix) {
  d <- cohort %>% mutate(trajectory = .data[[class_var]])
  names(d)[names(d) == "trajectory"] <- prefix
  fits <- list(
    list("Progression-free survival", "Unadjusted",
         as.formula(paste("Surv(pfs_time, pfs_event) ~", prefix))),
    list("Progression-free survival", "Adjusted",
         as.formula(paste("Surv(pfs_time, pfs_event) ~", prefix,
                          "+ age + sex + pathway"))),
    list("Overall survival", "Unadjusted",
         as.formula(paste("Surv(os_time, os_event) ~", prefix))),
    list("Overall survival", "Adjusted",
         as.formula(paste("Surv(os_time, os_event) ~", prefix,
                          "+ age + sex + pathway")))
  )
  bind_rows(lapply(fits, function(f) {
    s <- summary(coxph(f[[3]], data = d))
    keep <- startsWith(rownames(s$coefficients), prefix)
    tibble(outcome = f[[1]], model = f[[2]],
           term = rownames(s$coefficients)[keep],
           estimate = s$conf.int[keep, "exp(coef)"],
           conf.low = s$conf.int[keep, "lower .95"],
           conf.high = s$conf.int[keep, "upper .95"],
           p.value = s$coefficients[keep, "Pr(>|z|)"])
  }))
}

dashboard_outputs <- list(
  full_cohort = list(summary = tibble(
    patients = n, median_age = median(cohort$age),
    q1_age = quantile(cohort$age, 0.25), q3_age = quantile(cohort$age, 0.75),
    r2iss_available = sum(cohort$r2iss != "Missing")
  )),
  pp  = list(class_summary = class_summary("class_a", "trajectory_class"),
             cox_results   = cox_results("class_a", a_prefix)),
  b2m = list(class_summary = class_summary("class_b", "trajectory_b2m_class"),
             cox_results   = cox_results("class_b", b_prefix))
)
saveRDS(dashboard_outputs, file.path(data_dir, "dashboard_outputs.rds"))


# 5. STATIC IMAGES ------------------------------------------------------------

save_png <- function(plot, name) {
  ggsave(file.path(image_dir, name), plot, width = 7, height = 5, dpi = 120,
         bg = "white")
}

bar_plot <- function(var, title) {
  ggplot(count(cohort, value = .data[[var]]), aes(value, n)) +
    geom_col(fill = "#1565C0", width = 0.6) +
    labs(title = title, subtitle = "Synthetic demo data", x = NULL,
         y = "Patients") +
    theme_minimal(base_size = 14)
}

save_png(ggplot(cohort, aes(age)) +
           geom_histogram(binwidth = 2, fill = "#1565C0", colour = "white") +
           labs(title = "Age at baseline", subtitle = "Synthetic demo data",
                x = "Age, years", y = "Patients") +
           theme_minimal(base_size = 14),
         "full_cohort_age.png")
save_png(bar_plot("sex", "Sex"), "full_cohort_sex.png")
save_png(bar_plot("r2iss", "Disease stage"), "full_cohort_stage.png")
save_png(bar_plot("pathway", "Treatment pathway"), "full_cohort_pathway.png")

km_plot <- function(time, event, class_var, labels, title) {
  fit <- survfit(as.formula(paste0("Surv(", time, ", ", event, ") ~ ",
                                   class_var)), data = cohort)
  km <- tibble(time = fit$time, surv = fit$surv,
               class = rep(labels, fit$strata))
  ggplot(km, aes(time, surv, colour = class)) +
    geom_step(linewidth = 1) +
    scale_colour_manual(values = setNames(blues[seq_along(labels)], labels)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = title, subtitle = "Synthetic demo data",
         x = "Months from landmark", y = "Survival probability",
         colour = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

save_png(km_plot("pfs_time", "pfs_event", "class_a", a_labels,
                 "Biomarker A: progression-free survival"),
         "biomarker_a_pfs.png")
save_png(km_plot("os_time", "os_event", "class_a", a_labels,
                 "Biomarker A: overall survival"),
         "biomarker_a_os.png")
save_png(km_plot("pfs_time", "pfs_event", "class_b", b_labels,
                 "Biomarker B: progression-free survival"),
         "biomarker_b_pfs.png")
save_png(km_plot("os_time", "os_event", "class_b", b_labels,
                 "Biomarker B: overall survival"),
         "biomarker_b_os.png")

message("Synthetic demo inputs written to dashboard/dashboard_data and ",
        "dashboard/www. Run shiny::runApp(\"dashboard\") to open the app.")
