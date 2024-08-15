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
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
affymap219 <- affymap219 %>% tibble()
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")


# DEFINE THE SET ####
set01 <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]


# DEFINE SOX9 PROBES ####
probes_sox9 <- affymap219 %>%
  dplyr::filter(Symb == "SOX9") %>%
  pull(AffyID)


# WRANGLE THE EXPRESSION DATA ####
df_exprs <- data_expressionset_k1208[
  featureNames(data_expressionset_k1208) %in% probes_sox9,
  data_expressionset_k1208$Patient %nin% c(15, 18)
] %>%
  exprs() %>%
  t() %>%
  as_tibble(rownames = "CEL") %>%
  rename_if(is.numeric, ~ paste("sox9_", ., sep = ""))


# WRANGLE THE PHENOTYPE DATA ####
df_pheno <- data_expressionset_k1208 %>%
  pData %>% 
  dplyr::select(CEL, Center, Patient, Felzartamab, Group) %>%
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


# MERGE THE DATA ####
df <- df_pheno %>%
  left_join(df_exprs, by = "CEL")



# FIT MODELS ####
# 11720447_s_at
fit_11720447_s_at <- lmer(sox9_11720447_s_at ~ Felzartamab * time + (time | Patient), df)
# summary(fit_11720447_s_at)
# performance::check_model(fit_11720447_s_at)
# check_model(fit_11720447_s_at)
# check_normality(fit_11720447_s_at)
# Get predictions
effect_plot_11720447_s_at <- ggeffects::ggpredict(fit_11720447_s_at, c("time", "Felzartamab"))
# 11720448_at
fit_11720448_at <- lmer(sox9_11720448_at ~ Felzartamab * time + (time | Patient), df)
# summary(fit_11720448_at)
# check_model(fit_11720448_at)
# check_normality(fit_11720448_at)
# Get predictions
effect_plot_11720448_at <- ggeffects::ggpredict(fit_11720448_at, c("time", "Felzartamab"))
# 11758216_s_at
fit_11758216_s_at <- lmer(sox9_11758216_s_at ~ Felzartamab * time + (time | Patient), df)
# summary(fit_11758216_s_at)
# check_model(fit_11758216_s_at)
# check_normality(fit_11758216_s_at)
# Get predictions
effect_plot_11758216_s_at <- ggeffects::ggpredict(fit_11758216_s_at, c("time", "Felzartamab"))


# EXTRACT THE SLOPES FROM THE MODEL FITS ####
# 11720447_s_at
# Create list object into which to place model results
slopes_11720447_s_at <- list()
# Acute slopes
slopes_11720447_s_at[[1]] <- get_slope(
  .model_obj = fit_11720447_s_at, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_11720447_s_at[[2]] <- get_slope(
  .model_obj = fit_11720447_s_at, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_11720447_s_at[[3]] <- get_slope(
  .model_obj = fit_11720447_s_at, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_11720447_s_at <- slopes_11720447_s_at %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")  %>% round(3)) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# 11720448_at
# Create list object into which to place model results
slopes_11720448_at <- list()
# Acute slopes
slopes_11720448_at[[1]] <- get_slope(
  .model_obj = fit_11720448_at, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_11720448_at[[2]] <- get_slope(
  .model_obj = fit_11720448_at, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_11720448_at[[3]] <- get_slope(
  .model_obj = fit_11720448_at, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_11720448_at <- slopes_11720448_at %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr") %>% round(3)) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# 11758216_s_at
# Create list object into which to place model results
slopes_11758216_s_at <- list()
# Acute slopes
slopes_11758216_s_at[[1]] <- get_slope(
  .model_obj = fit_11758216_s_at, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes_11758216_s_at[[2]] <- get_slope(
  .model_obj = fit_11758216_s_at, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes_11758216_s_at[[3]] <- get_slope(
  .model_obj = fit_11758216_s_at, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)
model_11758216_s_at <- slopes_11758216_s_at %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr") %>% round(3)) %>%
  dplyr::select(name, result, p_value, p_adjusted)


# CREATE TABLE SUMMARIES OF SLOPES ####
size_font <- 10
table_slope_11720447_s_at <- model_11720447_s_at %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "11720447_s_at Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_11720447_s_at %>% print(preview = "pptx")

table_slope_11720448_at <- model_11720448_at %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "11720448_at Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_11720447_s_at %>% print(preview = "pptx")

table_slope_11758216_s_at <- model_11758216_s_at %>%
  dplyr::select(-p_value) %>%
  flextable::qflextable() %>%
  flextable::set_table_properties(layout = "autofit") %>%
  flextable::set_header_labels(
    name = "Group", result = "11758216_s_at Slope (95%CI)",
    "p_adjusted" = "FDR"
  ) %>%
  flextable::bold(i = NULL, part = "header") %>%
  flextable::fontsize(size = size_font, part = "all")
# table_slope_11720447_s_at %>% print(preview = "pptx")


# WRANGLE RESULTS FOR PLOTTING ####
resuls_sox9 <- tibble(
  variable = c("11720447_s_at", "11720448_at", "11758216_s_at"),
  reference_data = df_pheno  %>% list,
  plot_data = list(effect_plot_11720447_s_at, effect_plot_11720448_at, effect_plot_11758216_s_at),
  slope_tables = list(table_slope_11720447_s_at, table_slope_11720448_at, table_slope_11758216_s_at)

)


# EXPORT RESULTS FOR PLOTTING ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(resuls_sox9, file = paste(saveDir, "results_sox9_k1208.RData", sep = ""))


# MAKE PLOTS OF 11720447_s_at MODEL RESULTS ####
plot_11720447_s_at <- effect_plot_11720447_s_at %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = df,
    aes(x = time, y = sox9_11720447_s_at, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = df,
    aes(x = time, y = sox9_11720447_s_at, color = factor(Felzartamab), group = Patient, linetype = linetype),
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
  scale_y_continuous(limits = c(4, 11)) +
  labs(
    # title = "D",
    x = "Time post-treatment (weeks)",
    y = "SOX9 expression\n(11720447_s_at) "
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


# MAKE PLOTS OF 11720448_at MODEL RESULTS ####
plot_11720448_at <- effect_plot_11720448_at %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = df,
    aes(x = time, y = sox9_11720448_at, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = df,
    aes(x = time, y = sox9_11720448_at, color = factor(Felzartamab), group = Patient, linetype = linetype),
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
    y = "SOX9 expression\n(11720448_at)"
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
  scale_y_continuous(limits = c(4, 11)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )


# MAKE PLOTS OF 11758216_s_at MODEL RESULTS ####
plot_11758216_s_at <- effect_plot_11758216_s_at %>%
  ggplot(aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  geom_point(
    data = df,
    aes(x = time, y = sox9_11758216_s_at, color = factor(Felzartamab), group = Patient),
    alpha = 0,
    size = 0.1
  ) +
  geom_line(
    data = df,
    aes(x = time, y = sox9_11758216_s_at, color = factor(Felzartamab), group = Patient, linetype = linetype),
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
    y = "SOX9 expression\n(11758216_s_at)"
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
  scale_y_continuous(limits = c(4, 11)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )



# GENERATE GROBS FOR SLOPE TABLES ####
table_slope_11720447_s_at_grob <- table_slope_11720447_s_at %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()
table_slope_11720448_at_grob <- table_slope_11720448_at %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()
table_slope_11758216_s_at_grob <- table_slope_11758216_s_at %>%
  gen_grob(fit = "fixed", just = "centre") %>%
  as_ggplot()



# EXTRACT LEGEND FOR REGRESSION PLOTS ####
panel_legend_regression <- plot_11720447_s_at %>%
  ggpubr::get_legend() %>%
  ggpubr::as_ggplot()


# MAKE PANEL OF SLOPE PLOTS ####
plot_panels_regression <- patchwork::wrap_plots(
  plot_11720447_s_at, plot_11720448_at, plot_11758216_s_at,
  table_slope_11720447_s_at_grob, table_slope_11720448_at_grob, table_slope_11758216_s_at_grob,
  nrow = 2,
  ncol = 3,
  guides = "collect",
  heights = c(1, 0.5)
) +
  patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", rep("", 3)))) &
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    plot.tag = element_text(size = 20, face = "bold", vjust = 1)
  )

plot_panels_SOX9 <- ggarrange(
  panel_legend_regression,
  plot_panels_regression,
  nrow = 2,
  heights = c(0.15, 1.5)
)



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
  filename = paste(saveDir, "Felzartamab SOX9 regression 1208.png"),
  plot = plot_panels_SOX9,
  dpi = 300,
  width = 40,
  height = 15,
  units = "cm",
  bg = "white"
)


