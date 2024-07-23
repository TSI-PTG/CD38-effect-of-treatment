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
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/get_slope_function.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")
# load violin plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_plots.RData")
# load volcano enrichment plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano enrichment plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")


# WRANGLE THE DATA ####
m <- data_scores_k1208 %>%
  dplyr::select(Center, Patient, Felzartamab, Group, IRRAT30, IRITD3, IRITD5) %>%
  mutate(
    time = case_when(
      Group == "Index" ~ 0,
      Group == "FU1" ~ 0.46,
      Group == "FU2" ~ 0.99
    ),
    # linetype = ifelse(Patient == 9, "solid", "dashed") %>% factor()
    linetype = "solid" %>% factor()
  ) %>%
  dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18))



# FIT MODELS ####
# IRRAT30
fit_IRRAT30 <- lmer(IRRAT30 ~ Felzartamab * time + (time | Patient), m)
# summary(fit_IRRAT30)
# performance::check_model(fit_IRRAT30)
# check_model(fit_IRRAT30)
# check_normality(fit_IRRAT30)
# Get predictions
effect_plot_IRRAT30 <- ggeffects::ggpredict(fit_IRRAT30, c("time", "Felzartamab"))
# IRITD3
fit_IRITD3 <- lmer(IRITD3 ~ Felzartamab * time + (time | Patient), m)
# summary(fit_IRITD3)
# check_model(fit_IRITD3)
# check_normality(fit_IRITD3)
# Get predictions
effect_plot_IRITD3 <- ggeffects::ggpredict(fit_IRITD3, c("time", "Felzartamab"))
# IRITD5
fit_IRITD5 <- lmer(IRITD5 ~ Felzartamab * time + (time | Patient), m)
# summary(fit_IRITD5)
# check_model(fit_IRITD5)
# check_normality(fit_IRITD5)
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
size_font <- 10
table_slope_IRRAT30 <- model_IRRAT30 %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "IRRAT30 Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_IRRAT30 %>% print(preview = "pptx")

table_slope_IRITD3 <- model_IRITD3 %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "IRITD3 Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_IRRAT30 %>% print(preview = "pptx")

table_slope_IRITD5 <- model_IRITD5 %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "IRITD5 Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_IRRAT30 %>% print(preview = "pptx")


# WRANGLE RESULTS FOR PLOTTING ####
resuls_injury <- tibble(
  variable = c("IRRAT30", "IRITD3", "IRITD5"),
  reference_data = m  %>% list,
  plot_data = list(effect_plot_IRRAT30, effect_plot_IRITD3, effect_plot_IRITD5),
  slope_tables = list(table_slope_IRRAT30, table_slope_IRITD3, table_slope_IRITD5)

)

# EXPORT RESULTS FOR PLOTTING ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(resuls_injury, file = paste(saveDir, "results_injury_k1208.RData", sep = ""))


# MAKE PLOTS OF IRRAT30 MODEL RESULTS ####
plot_IRRAT30 <- effect_plot_IRRAT30 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = m,
    aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = m,
    aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient, linetype = linetype),
    stat = "smooth",
    method = "lm",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  # annotate(
  #   geom = "text",
  #   label = "Slope difference: -0.60 (95%CI -1.10 to -0.10), p=0.055",
  #   x = 0.5, y = 1.8, size = 4
  # ) +
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
    # title = "D",
    x = "Time post-treatment (weeks)",
    y = "Injury-repair associated\n(IRRAT30)"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.width = unit(4, "cm"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )


# MAKE PLOTS OF IRITD3 MODEL RESULTS ####
plot_IRITD3 <- effect_plot_IRITD3 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = m,
    aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = m,
    aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient, linetype = linetype),
    stat = "smooth",
    method = "lm",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  # annotate(
  #   geom = "text",
  #   label = "Slope difference: -0.20 (95%CI -0.37 to -0.02), p=0.041",
  #   x = 0.5, y = 0.65, size = 4
  # ) +
  labs(
    # title = "E",
    x = "Time post-treatment (weeks)",
    y = "Injury-repair induced, day 3\n(IRITD3)"
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
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )


# MAKE PLOTS OF IRITD5 MODEL RESULTS ####
plot_IRITD5 <- effect_plot_IRITD5 %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = m,
    aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = m,
    aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient, linetype = linetype),
    stat = "smooth",
    method = "lm",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  # annotate(
  #   geom = "text",
  #   label = "Slope difference: -0.23 (95%CI -0.41 to -0.05), p=0.042",
  #   x = 0.5, y = 0.75, size = 4
  # ) +
  labs(
    # title = "F",
    x = "Time post-treatment (weeks)",
    y = "Injury-repair induced, day 5\n(IRITD5)"
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
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )


# COMBINE THE SLOPE TABLES ####
# injury_slopes <- rbind(model_IRRAT30, model_IRITD3, model_IRITD5) %>%
#   tibble::add_row(name = "IRRAT", .after = 0) %>%
#   tibble::add_row(name = "IRITD3", .after = 4) %>%
#   tibble::add_row(name = "IRITD5", .after = 8)
# injury_slopes <- qflextable(injury_slopes) %>%
#   set_table_properties(layout = "autofit") %>%
#   set_header_labels(
#     name = "Scores", result = "Estimated slope (95%CI)",
#     "p_value" = "FDR", "p_adjusted" = "Adjusted FDR"
#   ) %>%
#   bold(i = NULL, part = "header") %>%
#   bold(i = c(1, 5, 9), j = 1) %>%
#   padding(i = c(2:4, 6:8, 10:12), j = 1, padding.left = 20)
# Print tables
# injury_slopes %>% print(preview = "pptx")


# GENERATE GROBS FOR SLOPE TABLES ####
table_slope_IRRAT30_grob <- table_slope_IRRAT30 %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()
table_slope_IRITD3_grob <- table_slope_IRITD3 %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()
table_slope_IRITD5_grob <- table_slope_IRITD5 %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()


# EXTRACT LEGEND FOR VIOLIN PLOTS ####
panel_legend_violin <- felzartamab_plots %>%
  dplyr::filter(variable == "AMAT1") %>%
  pull(plot_violin) %>%
  ggpubr::get_legend() %>%
  ggpubr::as_ggplot() +
  theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# EXTRACT LEGEND FOR REGRESSION PLOTS ####
panel_legend_regression <- plot_IRRAT30 %>%
  ggpubr::get_legend() %>%
  ggpubr::as_ggplot()


# EXTRACT JOINT LEGEND FOR REGRESSION AND VIOLIN PLOTS ####
panel_legend_all_injury <- panel_legend_violin + panel_legend_regression



# EXTRACT INJURY VIOLIN PLOTS ####
plots_violin <- felzartamab_plots %>%
  dplyr::filter(category %in% c("injury")) %>%
  pull(plot_violin)


# MAKE PANEL OF SLOPE PLOTS ####
plot_panels_regression <- patchwork::wrap_plots(
  plot_IRRAT30, plot_IRITD3, plot_IRITD5,
  table_slope_IRRAT30_grob, table_slope_IRITD3_grob, table_slope_IRITD5_grob,
  nrow = 2,
  ncol = 3,
  guides = "collect",
  heights = c(1, 0.5)
) +
  patchwork::plot_annotation(tag_levels = list(c("D", "E", "F", rep("", 3)))) &
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    plot.tag = element_text(size = 20, face = "bold", vjust = 1)
  )

plot_panels_regression_legend <- ggarrange(
  panel_legend_regression,
  plot_panels_regression,
  nrow = 2,
  heights = c(0.15, 1.5)
)


# MAKE PANEL OF VIOLIN PLOTS ####
# plot_panel_violin <- plots_violin %>%
#   wrap_plots(nrow = 1, ncol = 3) +
#   plot_annotation(tag_levels = list(c(LETTERS[1:15]))) &
#   theme(
#     legend.position = "none",
#     axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
#   )

# panels_violin_legend <- ggarrange(
#   panel_legend_violin,
#   plot_panel_violin,
#   nrow = 2,
#   heights = c(0.25, 1.5)
# )


# MAKE PANEL OF SLOPE AND VIOLIN PLOTS ####
plot_panels_all_injury <- patchwork::wrap_plots(
  plots_violin[[1]], plots_violin[[2]], plots_violin[[3]],
  plot_IRRAT30, plot_IRITD3, plot_IRITD5,
  table_slope_IRRAT30_grob, table_slope_IRITD3_grob, table_slope_IRITD5_grob,
  nrow = 3,
  ncol = 3,
  guides = "collect",
  heights = c(1, 1, 0.5)
) +
  patchwork::plot_annotation(tag_levels = list(c(LETTERS[1:6], rep("", 3)))) &
  theme(
    legend.position = "none",
    # plot.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black"),
    plot.tag = element_text(size = 20, face = "bold", vjust = 1)
  )

panels_all_injury <- ggarrange(
  panel_legend_all_injury,
  plot_panels_all_injury,
  nrow = 2,
  heights = c(0.15, 1.5)
)


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
# ggsave(
#   filename = paste(saveDir, "Felzartamab injury regression.png"),
#   plot = plot_panels_regression_legend,
#   dpi = 300,
#   width = 40,
#   height = 15,
#   units = "cm",
#   bg = "white"
# )

# ggsave(
#   filename = paste(saveDir, "Felzartamab injury ART.png"),
#   plot = panels_violin_legend,
#   dpi = 300,
#   width = 36,
#   height = 11,
#   units = "cm",
#   bg = "white"
# )

ggsave(
  filename = paste(saveDir, "Felzartamab injury ART and regression.png"),
  plot = panels_all_injury,
  dpi = 300,
  width = 40,
  height = 26,
  units = "cm",
  bg = "white"
)


## Save plot and arrange with violin plots =====
# ggarrange(panel_pairs_injury,
# ggarrange(plot_IRRAT30, plot_IRITD3, plot_IRITD5,
# nrow = 1, common.legend = T), nrow = 2, ncol = 1)

### Save Tiff

# ggsave(here::here("graphics", "injury_slopes.tiff"),
#        width = 15, height = 9,
#        dpi = 300, compression = "lzw", type="cairo" )
