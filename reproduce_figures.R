#!/usr/bin/env Rscript

# This script includes the code to generate Figure1 to Figure5 of 'Who pays for the grid?' submitted to Nature Energy for review by Clark, Halsey, and Alvarez.

required_packages <- c(
  "haven", "dplyr", "tidyr", "ggplot2", "purrr", "broom", "nnet",
  "marginaleffects", "patchwork", "survey", "scales", "sandwich"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(broom)
  library(nnet)
  library(marginaleffects)
  library(patchwork)
  library(survey)
  library(scales)
  library(sandwich)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  getwd()
}
data_file <- file.path(script_dir, "nature_energy_replication_data.sav")
if (!file.exists(data_file)) {
  stop(
    "Missing ", data_file,
    ". Run create_replication_data.R once to create the reduced dataset."
  )
}

raw_data <- read_sav(data_file)

recode_baseline_outcomes <- function(data, include_q20 = FALSE) {
  data$Q16 <- relevel(
    factor(data$Q16, levels = c("1", "2", "3"),
           labels = c("Nat. Ban", "No Nat. Ban", "Not Sure N")),
    ref = "Not Sure N"
  )
  data$Q17 <- if_else(is.na(data$Q17), 1, data$Q17)
  data$Q17 <- relevel(
    factor(data$Q17, levels = c("1", "2", "3"),
           labels = c("State Ban", "No State Ban", "Not Sure S")),
    ref = "Not Sure S"
  )
  data$Q18 <- if_else(is.na(data$Q18), 1, data$Q18)
  data$Q18 <- relevel(
    factor(data$Q18, levels = c("1", "2", "3"),
           labels = c("County Ban", "No County Ban", "Not Sure C")),
    ref = "Not Sure C"
  )
  data$Q19 <- relevel(
    factor(data$Q19, levels = c("1", "2", "3", "4", "5"),
           labels = c("Support", "Support", "Oppose", "Oppose", "Not Sure L")),
    ref = "Not Sure L"
  )
  if (include_q20) {
    data$Q20 <- relevel(
      factor(data$Q20, levels = c("1", "2", "3", "4", "5"),
             labels = c("Water", "Electricity", "Noise", "Groundwater", "No concerns")),
      ref = "No concerns"
    )
  }
  data
}

# ---------------------------------------------------------------------------
# Figure 1: raw responses under the four initial treatments
# ---------------------------------------------------------------------------
message("Creating Figure1.pdf")

figure1_data <- recode_baseline_outcomes(raw_data, include_q20 = TRUE) %>%
  mutate(
    treat1 = as.integer(randomization_2 == 1),
    treat2 = as.integer(randomization_2 == 2),
    treat3 = as.integer(randomization_2 == 3),
    treat4 = as.integer(randomization_2 == 4)
  )

treatment_labels <- c(
  treat1 = "Control 1",
  treat2 = "Control 2",
  treat3 = "Context Usage 1\n(Gallons/GwH)",
  treat4 = "Context Usage 2\n(Percentage)"
)
treatment_order <- unname(treatment_labels)

get_pct_tables <- function(data, term, vars) {
  bind_rows(lapply(vars, function(v) {
    as.data.frame(table(data[[v]]) * 100 / nrow(data))
  })) %>%
    mutate(term = term)
}

get_unweighted_prop_ci <- function(data, responsevar) {
  data %>%
    filter(!is.na(.data[[responsevar]])) %>%
    group_by(response = .data[[responsevar]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      total = sum(count),
      prop_test = map2(count, total, ~ broom::tidy(prop.test(.x, .y)))
    ) %>%
    unnest_wider(prop_test) %>%
    select(response, estimate, conf.low, conf.high) %>%
    rename(Var1 = response)
}

figure1_vars <- c("Q17", "Q16", "Q18", "Q19")
figure1_main <- bind_rows(lapply(names(treatment_labels), function(tr) {
  get_pct_tables(
    figure1_data %>% filter(.data[[tr]] == 1) %>% select(all_of(figure1_vars)),
    tr,
    figure1_vars
  )
})) %>%
  mutate(
    intervention = treatment_labels[term],
    experiment = case_when(
      Var1 %in% c("Nat. Ban", "No Nat. Ban", "Not Sure N") ~ "National Ban",
      Var1 %in% c("State Ban", "No State Ban", "Not Sure S") ~ "State Ban",
      Var1 %in% c("County Ban", "No County Ban", "Not Sure C") ~ "County Ban",
      TRUE ~ "Local Datacenter"
    ),
    stance = case_when(
      Var1 %in% c("No Nat. Ban", "No State Ban", "No County Ban", "Support") ~ "Pro-",
      Var1 %in% c("Nat. Ban", "State Ban", "County Ban", "Oppose") ~ "Anti-",
      TRUE ~ "Ambivalent"
    )
  )

figure1_ci <- bind_rows(lapply(names(treatment_labels), function(tr) {
  bind_rows(lapply(figure1_vars, function(v) {
    get_unweighted_prop_ci(filter(figure1_data, .data[[tr]] == 1), v) %>%
      mutate(term = tr)
  }))
}))

figure1_main <- figure1_main %>%
  left_join(figure1_ci, by = c("Var1", "term")) %>%
  mutate(
    conf.low = 100 * conf.low,
    conf.high = 100 * conf.high,
    experiment = factor(
      experiment,
      levels = c("National Ban", "State Ban", "County Ban", "Local Datacenter")
    ),
    intervention = factor(intervention, levels = treatment_order),
    Var1 = factor(
      Var1,
      levels = c(
        "Not Sure N", "No Nat. Ban", "Nat. Ban",
        "Not Sure S", "No State Ban", "State Ban",
        "Not Sure C", "No County Ban", "County Ban",
        "Not Sure L", "Support", "Oppose"
      )
    )
  )

figure1_q20 <- bind_rows(lapply(names(treatment_labels), function(tr) {
  get_pct_tables(
    figure1_data %>% filter(.data[[tr]] == 1) %>% select(Q20),
    tr,
    "Q20"
  )
})) %>%
  mutate(intervention = factor(treatment_labels[term], levels = treatment_order))

figure1_q20_ci <- bind_rows(lapply(names(treatment_labels), function(tr) {
  get_unweighted_prop_ci(filter(figure1_data, .data[[tr]] == 1), "Q20") %>%
    mutate(term = tr)
}))

figure1_q20 <- figure1_q20 %>%
  left_join(figure1_q20_ci, by = c("Var1", "term")) %>%
  mutate(
    conf.low = 100 * conf.low,
    conf.high = 100 * conf.high,
    Var1 = factor(
      Var1,
      levels = c("No concerns", "Noise", "Groundwater", "Water", "Electricity")
    )
  )

figure1_panel_a <- ggplot(figure1_main, aes(fill = stance, y = Var1, x = Freq)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.15, color = "grey40", orientation = "y"
  ) +
  geom_text(
    aes(label = sprintf("%.1f", Freq)),
    size = 2, position = position_nudge(y = 0.28), show.legend = FALSE
  ) +
  labs(title = "Panel A", y = "", x = "Percentage", fill = "") +
  scale_fill_manual(values = c(Ambivalent = "grey", `Pro-` = "red", `Anti-` = "blue")) +
  geom_vline(xintercept = 50, lty = 2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 270, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 12)
  ) +
  facet_grid(rows = vars(experiment), cols = vars(intervention), scales = "free_y")

figure1_panel_b <- ggplot(figure1_q20, aes(fill = Var1, y = Var1, x = Freq)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.15, color = "grey40", orientation = "y"
  ) +
  geom_text(
    aes(label = sprintf("%.1f", Freq)),
    size = 2, position = position_nudge(y = 0.28), show.legend = FALSE
  ) +
  labs(title = "Panel B", y = "", x = "Percentage", fill = "") +
  scale_fill_manual(values = c(
    Water = "skyblue", Electricity = "orange", Noise = "plum",
    Groundwater = "brown", `No concerns` = "forestgreen"
  )) +
  geom_vline(xintercept = 50, lty = 2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 270, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 12)
  ) +
  facet_grid(cols = vars(intervention))

figure1 <- figure1_panel_a / figure1_panel_b + plot_layout(heights = c(2.2, 1))
ggsave(file.path(script_dir, "Figure1.pdf"), figure1, width = 8, height = 11, units = "in")

# ---------------------------------------------------------------------------
# Figure 2: baseline demographic regressions
# ---------------------------------------------------------------------------
message("Creating Figure2.pdf")

figure2_data <- recode_baseline_outcomes(raw_data) %>%
  mutate(
    age4 = relevel(
      factor(age4, levels = c("1", "2", "3", "4"),
             labels = c("Under 30", "30-44", "45-64", "65+")),
      ref = "Under 30"
    ),
    gender4 = relevel(
      factor(gender_caltech, levels = c("1", "2", "3", "4"),
             labels = c("Man", "Woman", "Other", "Other")),
      ref = "Man"
    ),
    pid3 = relevel(
      factor(pid3, levels = c("1", "2", "3", "4", "5"),
             labels = c("Democrat", "Republican", "Independent", "Other", "Other")),
      ref = "Republican"
    ),
    educ4 = relevel(
      factor(educ4, levels = c("1", "2", "3", "4"),
             labels = c("HS or less", "College (2yr)", "College (4yr)", "Post-grad")),
      ref = "HS or less"
    ),
    aiuse = relevel(
      factor(Q9, levels = c("1", "2", "3", "4", "5"),
             labels = c("High", "High", "High", "Rare", "Rare")),
      ref = "Rare"
    ),
    drilling = factor(Q2_1, levels = c("1", "2", "3"),
                      labels = c("Pro", "Anti", "Not Sure")),
    nuclear = factor(Q2_2, levels = c("1", "2", "3"),
                     labels = c("Pro", "Anti", "Not Sure")),
    coal = factor(Q2_3, levels = c("1", "2", "3"),
                  labels = c("Pro", "Anti", "Not Sure")),
    solar = factor(Q2_4, levels = c("1", "2", "3"),
                   labels = c("Pro", "Anti", "Not Sure")),
    fracking = factor(Q2_5, levels = c("1", "2", "3"),
                      labels = c("Pro", "Anti", "Not Sure")),
    wind = factor(Q2_6, levels = c("1", "2", "3"),
                  labels = c("Pro", "Anti", "Not Sure")),
    n_energy_supported =
      as.integer(drilling == "Pro") + as.integer(nuclear == "Pro") +
      as.integer(coal == "Pro") + as.integer(solar == "Pro") +
      as.integer(fracking == "Pro") + as.integer(wind == "Pro"),
    proenergy = relevel(
      factor(
        case_when(
          n_energy_supported < 3 ~ "Fewer than 3",
          n_energy_supported == 3 ~ "Exactly 3",
          n_energy_supported > 3 ~ "More than 3"
        ),
        levels = c("Fewer than 3", "Exactly 3", "More than 3")
      ),
      ref = "Fewer than 3"
    )
  )

extract_baseline_effects <- function(outcome) {
  model <- multinom(
    as.formula(paste0(
      outcome, " ~ age4 + gender4 + pid3 + educ4 + aiuse + proenergy"
    )),
    data = figure2_data,
    trace = FALSE
  )
  result <- as.data.frame(avg_comparisons(model, wts = model$weights))
  result$question <- outcome
  result
}

figure2_effects <- bind_rows(lapply(c("Q16", "Q17", "Q18", "Q19"), extract_baseline_effects)) %>%
  separate_wider_delim(contrast, " - ", names = c("comp", "ref")) %>%
  mutate(
    comp = gsub("^mean\\(|\\)$", "", comp),
    ref = gsub("^mean\\(|\\)$", "", ref)
  ) %>%
  filter(
    (question == "Q16" & group == "No Nat. Ban") |
      (question == "Q17" & group == "No State Ban") |
      (question == "Q18" & group == "No County Ban") |
      (question == "Q19" & group == "Support")
  ) %>%
  filter(!(comp %in% c("NA", "Other"))) %>%
  mutate(
    signif = factor(
      case_when(conf.low > 0 ~ "pos", conf.high < 0 ~ "neg", TRUE ~ "ns"),
      levels = c("pos", "ns", "neg")
    ),
    question = factor(
      case_when(
        question == "Q16" ~ "Against\nNational Ban",
        question == "Q17" ~ "Against\nState Ban",
        question == "Q18" ~ "Against\nCounty Ban",
        TRUE ~ "For\nData Center"
      ),
      levels = c(
        "Against\nNational Ban", "Against\nState Ban",
        "Against\nCounty Ban", "For\nData Center"
      )
    ),
    term_f = factor(
      case_when(
        term == "age4" ~ "Age\n(Under-30)",
        term == "gender4" ~ "Gender\n(Man)",
        term == "pid3" ~ "Party ID\n(Republican)",
        term == "educ4" ~ "Education\n(HS or less)",
        term == "aiuse" ~ "AI use\n(Rare)",
        term == "proenergy" ~ "Energy support\n(Fewer than 3)",
        TRUE ~ term
      ),
      levels = c(
        "Age\n(Under-30)", "Gender\n(Man)", "Party ID\n(Republican)",
        "Education\n(HS or less)", "AI use\n(Rare)",
        "Energy support\n(Fewer than 3)"
      )
    )
  )

figure2_min <- min(figure2_effects$conf.low, na.rm = TRUE)
figure2_max <- max(figure2_effects$conf.high, na.rm = TRUE)
figure2_pad <- 0.02 * (figure2_max - figure2_min)

figure2 <- ggplot(figure2_effects, aes(x = estimate, y = comp, colour = signif)) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.12, orientation = "y"
  ) +
  coord_cartesian(xlim = c(figure2_min - figure2_pad, figure2_max + figure2_pad)) +
  scale_color_manual(values = c(pos = "blue", ns = "grey60", neg = "red")) +
  labs(
    title = "Relative anti-ban and pro-data-center preference based on demographics",
    x = "Average marginal effect",
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  facet_grid(rows = vars(term_f), cols = vars(question), scales = "free_y", space = "free_y")

ggsave(file.path(script_dir, "Figure2.pdf"), figure2, width = 10.5, height = 13.5, units = "in")

# ---------------------------------------------------------------------------
# Figure 3: weighted treatment responses and effects
# ---------------------------------------------------------------------------
message("Creating Figure3.pdf")

recode_followup_outcomes <- function(data) {
  data$Q21 <- relevel(
    factor(data$Q21, levels = c("1", "2", "3"),
           labels = c("Nat. Ban", "No Nat. Ban", "Not Sure N")),
    ref = "Not Sure N"
  )
  data$Q22 <- if_else(is.na(data$Q22), 1, data$Q22)
  data$Q22 <- relevel(
    factor(data$Q22, levels = c("1", "2", "3"),
           labels = c("State Ban", "No State Ban", "Not Sure S")),
    ref = "Not Sure S"
  )
  data$Q23 <- if_else(is.na(data$Q23), 1, data$Q23)
  data$Q23 <- relevel(
    factor(data$Q23, levels = c("1", "2", "3"),
           labels = c("County Ban", "No County Ban", "Not Sure C")),
    ref = "Not Sure C"
  )
  data$Q24 <- relevel(
    factor(data$Q24, levels = c("1", "2", "3", "4", "5"),
           labels = c("Support", "Support", "Oppose", "Oppose", "Not Sure L")),
    ref = "Not Sure L"
  )
  data
}

figure3_data <- raw_data %>%
  mutate(
    randomization_3 = zap_labels(randomization_3),
    randomization_3 = case_when(
      randomization_3 == 2 ~ 3,
      randomization_3 == 3 ~ 2,
      TRUE ~ randomization_3
    )
  ) %>%
  recode_baseline_outcomes() %>%
  recode_followup_outcomes() %>%
  mutate(
    treat1 = as.integer(randomization_2 == 1),
    treat2 = as.integer(randomization_2 == 2),
    treat3 = as.integer(randomization_2 == 3),
    treat6 = as.integer(randomization_3 == 1),
    treat7 = as.integer(randomization_3 == 2),
    treat8 = as.integer(randomization_3 == 3)
  )

weighted_prop_ci <- function(data, responsevar) {
  analysis_data <- data %>%
    transmute(resp = .data[[responsevar]], w = as.numeric(weight)) %>%
    filter(!is.na(resp), !is.na(w), is.finite(w), w > 0)
  design <- svydesign(ids = ~1, weights = ~w, data = analysis_data)
  estimate <- svymean(~resp, design = design, na.rm = TRUE)
  interval <- suppressMessages(confint(estimate))
  tibble(
    Var1 = gsub("^resp", "", names(coef(estimate))),
    estimate = as.numeric(coef(estimate)),
    conf.low = interval[, 1],
    conf.high = interval[, 2]
  )
}

ban_relabel <- function(x) {
  case_when(
    x == "Nat. Ban" ~ "For Nat. Ban",
    x == "No Nat. Ban" ~ "Against Nat. Ban",
    x == "State Ban" ~ "For State Ban",
    x == "No State Ban" ~ "Against State Ban",
    x == "County Ban" ~ "For County Ban",
    x == "No County Ban" ~ "Against County Ban",
    TRUE ~ x
  )
}
not_sure_relabel <- function(x) {
  if_else(
    x %in% c("Not Sure N", "Not Sure S", "Not Sure C", "Not Sure L"),
    "Not Sure",
    x
  )
}

figure3_raw_parts <- list(
  list(data = figure3_data, questions = c("Q17", "Q16", "Q18", "Q19"), term = "treat1"),
  list(data = filter(figure3_data, treat6 == 1), questions = c("Q22", "Q21", "Q23", "Q24"), term = "treat2"),
  list(data = filter(figure3_data, treat7 == 1), questions = c("Q22", "Q21", "Q23", "Q24"), term = "treat3"),
  list(data = filter(figure3_data, treat8 == 1), questions = c("Q22", "Q21", "Q23", "Q24"), term = "treat4")
)

figure3_raw <- bind_rows(lapply(figure3_raw_parts, function(part) {
  bind_rows(lapply(part$questions, function(question) {
    weighted_prop_ci(part$data, question) %>% mutate(term = part$term)
  }))
})) %>%
  mutate(
    Freq = 100 * estimate,
    conf.low = 100 * conf.low,
    conf.high = 100 * conf.high,
    intervention = case_when(
      term == "treat1" ~ "(1)\nCombined\n Pre-treatment",
      term == "treat2" ~ "(2)\nProperty Tax",
      term == "treat3" ~ "(3)\nProperty Tax \n VA elec",
      TRUE ~ "(4)\nProperty Tax \n VA + ME elec"
    ),
    experiment = factor(
      case_when(
        Var1 %in% c("Nat. Ban", "No Nat. Ban", "Not Sure N") ~ "National Ban",
        Var1 %in% c("State Ban", "No State Ban", "Not Sure S") ~ "State Ban",
        Var1 %in% c("County Ban", "No County Ban", "Not Sure C") ~ "County Ban",
        TRUE ~ "Local Datacenter"
      ),
      levels = c("National Ban", "State Ban", "County Ban", "Local Datacenter")
    ),
    stance = case_when(
      Var1 %in% c("No Nat. Ban", "No State Ban", "No County Ban", "Support") ~ "Pro-",
      Var1 %in% c("Nat. Ban", "State Ban", "County Ban", "Oppose") ~ "Anti-",
      TRUE ~ "Ambivalent"
    ),
    Var1_plot = not_sure_relabel(ban_relabel(Var1))
  )

figure3_levels <- not_sure_relabel(ban_relabel(c(
  "Not Sure N", "No Nat. Ban", "Nat. Ban",
  "Not Sure S", "No State Ban", "State Ban",
  "Not Sure C", "No County Ban", "County Ban",
  "Not Sure L", "Support", "Oppose"
)))
figure3_raw$Var1_plot <- factor(
  figure3_raw$Var1_plot,
  levels = figure3_levels[!duplicated(figure3_levels)]
)

figure3_panel_a <- ggplot(figure3_raw, aes(fill = stance, y = Var1_plot, x = Freq)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.15, color = "grey40", orientation = "y"
  ) +
  geom_text(
    aes(label = sprintf("%.1f", Freq)),
    size = 2.1, position = position_nudge(y = 0.28), show.legend = FALSE
  ) +
  labs(title = "Responses under treatments", y = "", x = "Percentage", fill = "") +
  scale_fill_manual(values = c(Ambivalent = "grey", `Pro-` = "red", `Anti-` = "blue")) +
  geom_vline(xintercept = 50, lty = 2) +
  geom_vline(xintercept = 25, lty = 3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 270, size = 11),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 10.5),
    plot.title = element_text(face = "bold")
  ) +
  facet_grid(rows = vars(experiment), cols = vars(intervention), scales = "free_y")

figure3_controls <- figure3_data %>%
  mutate(treat6 = 0L, treat7 = 0L, treat8 = 0L) %>%
  select(-Q21, -Q22, -Q23, -Q24) %>%
  rename(Q21 = Q16, Q22 = Q17, Q23 = Q18, Q24 = Q19)
figure3_stacked <- bind_rows(figure3_data, figure3_controls)

fit_followup_effect <- function(outcome, treatment, randomization_value) {
  model <- multinom(
    as.formula(paste0(outcome, " ~ ", treatment, " + treat1 + treat2 + treat3")),
    data = figure3_stacked %>% filter(randomization_3 == randomization_value),
    trace = FALSE
  )
  avg_comparisons(model, variables = treatment, wts = model$weights)
}

figure3_effects <- bind_rows(lapply(c("Q21", "Q22", "Q23", "Q24"), function(outcome) {
  bind_rows(
    fit_followup_effect(outcome, "treat6", 1),
    fit_followup_effect(outcome, "treat7", 2),
    fit_followup_effect(outcome, "treat8", 3)
  )
})) %>%
  as.data.frame() %>%
  mutate(
    intervention = case_when(
      term == "treat6" ~ "(1)\nProperty Tax",
      term == "treat7" ~ "(2)\nProperty Tax \n VA elec",
      TRUE ~ "(3)\nProperty Tax \n VA elec + ME elec"
    ),
    experiment = factor(
      case_when(
        group %in% c("Nat. Ban", "No Nat. Ban", "Not Sure N") ~ "National Ban",
        group %in% c("State Ban", "No State Ban", "Not Sure S") ~ "State Ban",
        group %in% c("County Ban", "No County Ban", "Not Sure C") ~ "County Ban",
        TRUE ~ "Local Datacenter"
      ),
      levels = c("National Ban", "State Ban", "County Ban", "Local Datacenter")
    ),
    stance = case_when(
      group %in% c("No Nat. Ban", "No State Ban", "No County Ban", "Support") ~ "Pro-",
      group %in% c("Nat. Ban", "State Ban", "County Ban", "Oppose") ~ "Anti-",
      TRUE ~ "Ambivalent"
    ),
    signif = case_when(conf.high < 0 ~ "neg", conf.low > 0 ~ "pos", TRUE ~ "nosig"),
    stance = if_else(signif == "nosig", "Ambivalent", stance),
    stance = if_else(signif == "neg" & stance == "Pro-", "Anti-", stance),
    stance = if_else(signif == "neg" & stance == "Anti-", "Pro-", stance),
    group_plot = not_sure_relabel(ban_relabel(group))
  )
figure3_effects$group_plot <- factor(
  figure3_effects$group_plot,
  levels = figure3_levels[!duplicated(figure3_levels)]
)

figure3_panel_b <- ggplot(
  figure3_effects,
  aes(x = estimate, y = group_plot, colour = stance)
) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  coord_cartesian(xlim = range(c(
    figure3_effects$conf.low,
    figure3_effects$conf.high
  ), na.rm = TRUE)) +
  labs(
    title = "Treatment effects under treatments",
    x = "Average marginal effect",
    y = ""
  ) +
  scale_color_manual(values = c(Ambivalent = "grey40", `Anti-` = "red", `Pro-` = "blue")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 10.5),
    plot.title = element_text(face = "bold")
  ) +
  geom_errorbar(width = 0.1, aes(xmin = conf.low, xmax = conf.high)) +
  facet_grid(rows = vars(experiment), cols = vars(intervention), scales = "free_y")

figure3 <- (figure3_panel_a / figure3_panel_b) +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(file.path(script_dir, "Figure3.pdf"), figure3, width = 8, height = 11, units = "in")

# ---------------------------------------------------------------------------
# Shared conjoint data for Figures 4 and 5
# ---------------------------------------------------------------------------
conjoint_data <- raw_data %>%
  mutate(
    respondent_id = row_number(),
    prefer12 = as.integer(prefer12),
    prefer34 = as.integer(prefer34),
    prefer56 = as.integer(prefer56)
  ) %>%
  filter(!is.na(prefer12), !is.na(prefer34), !is.na(prefer56)) %>%
  pivot_longer(
    cols = matches("^(company|location|energy|water|utility)_rand_[1-6]$"),
    names_to = c(".value", "profile"),
    names_pattern = "(.*)_rand_(\\d+)"
  ) %>%
  mutate(
    profile = as.integer(profile),
    company = relevel(
      factor(company, levels = c("1", "2"), labels = c("OpenAI", "Amazon")),
      ref = "Amazon"
    ),
    location = relevel(
      factor(location, levels = c("3", "4"), labels = c("5 miles", "50 miles")),
      ref = "5 miles"
    ),
    energy = relevel(
      factor(energy, levels = c("5", "6"), labels = c("Renewable", "Natural gas")),
      ref = "Natural gas"
    ),
    water = relevel(
      factor(water, levels = c("7", "8"),
             labels = c("500K gal/day", "2 mil gal/day")),
      ref = "2 mil gal/day"
    ),
    utility = relevel(
      factor(utility, levels = c("9", "10"), labels = c("+10%", "No increase")),
      ref = "+10%"
    ),
    party = relevel(
      factor(as.character(pid3), levels = c("1", "2", "3", "4", "5"),
             labels = c("Democrat", "Republican", "Independent", "Other", "Other")),
      ref = "Republican"
    ),
    aiuse = relevel(
      factor(as.character(Q9), levels = c("1", "2", "3", "4", "5"),
             labels = c("High", "High", "High", "Rare", "Rare")),
      ref = "Rare"
    ),
    prefer = case_when(
      profile == 1 & prefer12 == 1 ~ 1,
      profile == 2 & prefer12 == 2 ~ 1,
      profile %in% c(1, 2) ~ 0,
      profile == 3 & prefer34 == 1 ~ 1,
      profile == 4 & prefer34 == 2 ~ 1,
      profile %in% c(3, 4) ~ 0,
      profile == 5 & prefer56 == 1 ~ 1,
      profile == 6 & prefer56 == 2 ~ 1,
      TRUE ~ 0
    )
  )

conjoint_terms <- c("company", "location", "energy", "water", "utility")
conjoint_label_order <- c(
  "OpenAI (vs Amazon)",
  "50 miles (vs 5 miles)",
  "Renewable (vs Natural gas)",
  "500K gal/day (vs 2 mil gal/day)",
  "No increase (vs +10%)"
)

prepare_amce <- function(result) {
  as.data.frame(result) %>%
    separate_wider_delim(contrast, " - ", names = c("comp", "ref")) %>%
    mutate(
      comp = gsub("^mean\\(|\\)$", "", comp),
      ref = gsub("^mean\\(|\\)$", "", ref),
      label = case_when(
        term == "company" ~ paste0(comp, " (vs Amazon)"),
        term == "location" ~ paste0(comp, " (vs 5 miles)"),
        term == "energy" ~ paste0(comp, " (vs Natural gas)"),
        term == "water" ~ paste0(comp, " (vs 2 mil gal/day)"),
        term == "utility" ~ paste0(comp, " (vs +10%)")
      ),
      label = factor(label, levels = conjoint_label_order)
    )
}

pooled_amce_plot <- function(data) {
  ggplot(data, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      width = 0.15, color = "black", orientation = "y"
    ) +
    geom_point(shape = 1, size = 5, color = "black", stroke = 1) +
    labs(title = "Pooled AMCE", x = "Average marginal effect", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold")
    )
}

group_amce_plot <- function(data, group_var, colors, title) {
  dodge <- position_dodge(width = 0.55)
  ggplot(data, aes(x = estimate, y = label, colour = .data[[group_var]])) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      width = 0.15, orientation = "y", position = dodge
    ) +
    geom_point(shape = 1, size = 4, position = dodge, stroke = 1) +
    scale_color_manual(values = colors, name = NULL) +
    labs(title = title, x = "Average marginal effect", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold"),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.8),
        color = NA
      ),
      legend.text = element_text(size = 13)
    )
}

# ---------------------------------------------------------------------------
# Figure 4: pooled and subgroup conjoint AMCEs
# ---------------------------------------------------------------------------
message("Creating Figure4.pdf")

figure4_pooled_model <- glm(
  prefer ~ company + location + energy + water + utility,
  family = binomial(),
  data = conjoint_data,
  weights = weight
)
figure4_pooled <- prepare_amce(avg_comparisons(
  figure4_pooled_model,
  vcov = ~ respondent_id
))

figure4_party_model <- glm(
  prefer ~ party * (company + location + energy + water + utility),
  family = binomial(),
  data = conjoint_data,
  weights = weight
)
figure4_party <- prepare_amce(avg_comparisons(
  figure4_party_model,
  variables = conjoint_terms,
  by = "party",
  vcov = ~ respondent_id
)) %>%
  filter(party %in% c("Democrat", "Republican"))

figure4_ai_model <- glm(
  prefer ~ aiuse * (company + location + energy + water + utility),
  family = binomial(),
  data = conjoint_data,
  weights = weight
)
figure4_ai <- prepare_amce(avg_comparisons(
  figure4_ai_model,
  variables = conjoint_terms,
  by = "aiuse",
  vcov = ~ respondent_id
)) %>%
  filter(aiuse %in% c("High", "Rare"))

figure4 <- pooled_amce_plot(figure4_pooled) /
  (
    group_amce_plot(
      figure4_party, "party",
      c(Democrat = "blue", Republican = "red"),
      "by Party"
    ) |
      group_amce_plot(
        figure4_ai, "aiuse",
        c(High = "forestgreen", Rare = "grey50"),
        "by AI use"
      )
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file.path(script_dir, "Figure4.pdf"), figure4, width = 12.65, height = 8.7, units = "in")

# ---------------------------------------------------------------------------
# Figure 5: predictions for fixed conjoint profiles
# ---------------------------------------------------------------------------
message("Creating Figure5.pdf")

figure5_model <- figure4_pooled_model
figure5_vcov <- vcovCL(figure5_model, cluster = ~ respondent_id)

average_profile_prediction <- function(specification) {
  newdata <- conjoint_data
  for (variable in names(specification)) {
    newdata[[variable]] <- factor(
      specification[[variable]],
      levels = levels(conjoint_data[[variable]]),
      labels = levels(conjoint_data[[variable]])
    )
  }
  matrix <- model.matrix(delete.response(terms(figure5_model)), newdata)
  coefficients <- coef(figure5_model)
  prediction <- as.vector(plogis(matrix %*% coefficients))
  analysis_weights <- newdata$weight
  estimate <- weighted.mean(prediction, analysis_weights, na.rm = TRUE)
  gradient <- colSums(
    matrix * (prediction * (1 - prediction) * analysis_weights),
    na.rm = TRUE
  ) / sum(analysis_weights, na.rm = TRUE)
  standard_error <- sqrt(drop(
    t(gradient) %*% figure5_vcov %*% gradient
  ))
  tibble(
    estimate = estimate,
    std.error = standard_error,
    conf.low = estimate - qnorm(0.975) * standard_error,
    conf.high = estimate + qnorm(0.975) * standard_error
  )
}

fixed_profiles <- tibble(
  profile = c(
    "10%+ elec. rates",
    "2 million gallons of water",
    "10%+ elec. rates\n2 million gallons of water",
    "No impact on elec. rates",
    "10%+ elec. rates\n500,000 gallons of water",
    "No impact on elec. rates\n500,000 gallons of water",
    "No impact on elec. rates\n500,000 gallons of water\n50 miles from the house",
    "No impact on elec. rates\n500,000 gallons of water\n50 miles from the house\n100% renewable energy",
    "No impact on elec. rates\n500,000 gallons of water\n50 miles from the house\n100% renewable energy\nAmazon"
  ),
  specification = list(
    list(utility = "+10%"),
    list(water = "2 mil gal/day"),
    list(utility = "+10%", water = "2 mil gal/day"),
    list(utility = "No increase"),
    list(utility = "+10%", water = "500K gal/day"),
    list(utility = "No increase", water = "500K gal/day"),
    list(utility = "No increase", water = "500K gal/day", location = "50 miles"),
    list(
      utility = "No increase", water = "500K gal/day",
      location = "50 miles", energy = "Renewable"
    ),
    list(
      utility = "No increase", water = "500K gal/day",
      location = "50 miles", energy = "Renewable", company = "Amazon"
    )
  )
)

figure5_predictions <- fixed_profiles %>%
  rowwise() %>%
  mutate(prediction = list(average_profile_prediction(specification))) %>%
  unnest(prediction) %>%
  ungroup() %>%
  mutate(
    profile = factor(profile, levels = rev(profile)),
    estimate = pmin(pmax(estimate, 0), 1),
    conf.low = pmin(pmax(conf.low, 0), 1),
    conf.high = pmin(pmax(conf.high, 0), 1),
    label = percent(estimate, accuracy = 0.1)
  )

figure5 <- ggplot(figure5_predictions, aes(x = estimate, y = profile)) +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.15, color = "black", orientation = "y"
  ) +
  geom_point(shape = 16, size = 2.2, color = "black") +
  geom_text(aes(label = label), nudge_y = 0.22, size = 2.9, color = "black") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Average Pr(chosen)",
    y = "",
    title = "Predicted choice probability for fixed data center profiles"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(script_dir, "Figure5.pdf"), figure5, width = 7.8, height = 10, units = "in")

message("Finished. Wrote Figure1.pdf through Figure5.pdf to ", script_dir)
