# Linear Slope Model Injury =====

# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(lme4) # install.packages("lme4")
library(lmerTest) # install.packages("lmerTest")
library(performance) # install.packages("performance")
library(ggeffects) # install.packages("ggeffects")
library(flextable) # install.packages("flextable")
library(ggpubr) # install.packages("ggpubr")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/get_slope_function.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_felzartamab_k1208.RData")


# WRANGLE THE DATA ####
m <- data_felzartamab_k1208 %>%
  mutate(
    time = case_when(
      Group == "Index" ~ 0,
      Group == "FU1" ~ 0.46,
      Group == "FU2" ~ 0.99
    )
  ) %>%
  dplyr::filter(Group != "FU1b", Patient %nin% c(9, 15, 18))


# FIT MODELS ####
# IRRAT30
fit_IRRAT30 <- lmer(IRRAT30 ~ Felzartamab * time + (time | Patient), m)
# summary(fit_IRRAT30)
# performance::check_model(fit_IRRAT30)
check_model(fit_IRRAT30)
check_normality(fit_IRRAT30)
# Get predictions
effect_plot_IRRAT30 <- ggeffects::ggpredict(fit_IRRAT30, c("time", "Felzartamab"))
# IRITD3
fit_IRITD3 <- lmer(IRITD3 ~ Felzartamab * time + (time | Patient), m)
# summary(fit_IRITD3)
check_model(fit_IRITD3)
check_normality(fit_IRITD3)
# Get predictions
effect_plot_IRITD3 <- ggeffects::ggpredict(fit_IRITD3, c("time", "Felzartamab"))
# IRITD5
fit_IRITD5 <- lmer(IRITD5 ~ Felzartamab * time + (time | Patient), m)
summary(fit_IRITD5)
check_model(fit_IRITD5)
check_normality(fit_IRITD5)
# Get predictions
effect_plot_IRITD5 <- ggeffects::ggpredict(fit_IRITD5, c("time", "Felzartamab"))


# EXTRACT THE SLOPES FROM THE MODEL FITS ####
# IRRAT30
# Create list object into which to place model results
slopes_IRRAT30 <- list()
# Acute slopes
slopes_IRRAT30[[1]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_IRRAT30[[2]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_IRRAT30[[3]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_IRRAT30 <- slopes_IRRAT30 %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# IRITD3
# Create list object into which to place model results
slopes_IRITD3 <- list()
# Acute slopes
slopes_IRITD3[[1]] <- get_slope(
  .model_obj = fit_IRITD3, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_IRITD3[[2]] <- get_slope(
  .model_obj = fit_IRITD3, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_IRITD3[[3]] <- get_slope(
  .model_obj = fit_IRITD3, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_IRITD3 <- slopes_IRITD3 %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# IRITD5
# Create list object into which to place model results
slopes_IRITD5 <- list()
# Acute slopes
slopes_IRITD5[[1]] <- get_slope(
  .model_obj = fit_IRITD5, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_IRITD5[[2]] <- get_slope(
  .model_obj = fit_IRITD5, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_IRITD5[[3]] <- get_slope(
  .model_obj = fit_IRITD5, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_IRITD5 <- slopes_IRITD5 %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  dplyr::select(name, result, p_value, p_adjusted)


# CREATE TABLE SUMMARIES OF SLOPES ####
table_slope_IRRAT30 <- qflextable(model_IRRAT30) %>%
  set_table_properties(layout = "autofit") %>%
  set_header_labels(
    name = "Group", result = "IRRAT30 Slope (95%CI)",
    "p_value" = "P-Value"
  ) %>%
  bold(i = NULL, part = "header")
# table_slope_IRRAT30 %>% print(preview = "pptx")

table_slope_IRITD3 <- qflextable(model_IRITD3) %>%
  set_table_properties(layout = "autofit") %>%
  set_header_labels(
    name = "Group", result = "IRITD3 Slope (95%CI)",
    "p_value" = "P-Value"
  ) %>%
  bold(i = NULL, part = "header")
# table_slope_IRRAT30 %>% print(preview = "pptx")

table_slope_IRITD5 <- qflextable(model_IRITD5) %>%
  set_table_properties(layout = "autofit") %>%
  set_header_labels(
    name = "Group", result = "IRITD5 Slope (95%CI)",
    "p_value" = "P-Value"
  ) %>%
  bold(i = NULL, part = "header")
# table_slope_IRRAT30 %>% print(preview = "pptx")


# MAKE PLOTS OF IRRAT30 MODEL RESULTS ####
plot_IRRAT30 <- effect_plot_IRRAT30 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(data = m, aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.60 (95%CI -1.10 to -0.10), p=0.055",
    x = 0.5, y = 1.8, size = 4
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.46, 0.99),
    labels = c("0", "24", "52")
  ) +
  scale_y_continuous(limits = c(-1.5, 2)) +
  labs(
    title = "D",
    x = "Weeks post-treatment",
    y = "Injury-repair associated transcripts\n(IRRAT30)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )


# MAKE PLOTS OF IRITD3 MODEL RESULTS ####
plot_IRITD3 <- effect_plot_IRITD3 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(data = m, aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.20 (95%CI -0.37 to -0.02), p=0.041",
    x = 0.5, y = 0.65, size = 4
  ) +
  labs(
    title = "E",
    x = "Weeks post-treatment",
    y = "Injury-repair induced transcripts, day 3\n(IRITD3)"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.46, 1),
    labels = c("0", "24", "52")
  ) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )


# MAKE PLOTS OF IRITD5 MODEL RESULTS ####
plot_IRITD5 <- effect_plot_IRITD5 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(data = m, aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.23 (95%CI -0.41 to -0.05), p=0.042",
    x = 0.5, y = 0.75, size = 4
  ) +
  labs(
    title = "F",
    x = "Weeks post-treatment",
    y = "Injury-repair induced transcripts, day 5\n(IRITD5)"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.46, 0.99),
    labels = c("0", "24", "52")
  ) +
  scale_y_continuous(limits = c(-0.20, 0.8)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )


# COMBINE THE SLOPE TABLES ####
injury_slopes <- rbind(model_IRRAT30, model_IRITD3, model_IRITD5) %>%
  tibble::add_row(name = "IRRAT", .after = 0) %>%
  tibble::add_row(name = "IRITD3", .after = 4) %>%
  tibble::add_row(name = "IRITD5", .after = 8)
injury_slopes <- qflextable(injury_slopes) %>%
  set_table_properties(layout = "autofit") %>%
  set_header_labels(
    name = "Scores", result = "Estimated slope (95%CI)",
    "p_value" = "P-Value", "p_adjusted" = "Adjusted P-Value"
  ) %>%
  bold(i = NULL, part = "header") %>%
  bold(i = c(1, 5, 9), j = 1) %>%
  padding(i = c(2:4, 6:8, 10:12), j = 1, padding.left = 20)
# Print tables 
injury_slopes %>% print(preview = "pptx")


# COMBINE THE SLOPE PLOTS ####







## Save plot and arrange with violin plots =====
# ggarrange(panel_pairs_injury,
# ggarrange(plot_IRRAT30, plot_IRITD3, IRITD5_plot,
# nrow = 1, common.legend = T), nrow = 2, ncol = 1)

### Save Tiff

# ggsave(here::here("graphics", "injury_slopes.tiff"),
#        width = 15, height = 9,
#        dpi = 300, compression = "lzw", type="cairo" )
