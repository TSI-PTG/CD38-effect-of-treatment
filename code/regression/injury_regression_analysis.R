# Linear Slope Model Injury =====

# HOUSEKEEPING
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





## IRRAT30 ====
fit_IRRAT30 <- lmer(IRRAT30 ~ Felzartamab * time + (time | Patient), m)
summary(fit_IRRAT30)
performance::check_model(fit_IRRAT30)

# Get predictions
effect_plot <- ggeffects::ggpredict(fit_IRRAT30, c("time", "Felzartamab"))

# Plot slopes
irrat_plot <-
  ggplot(effect_plot, aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  theme_classic() +
  labs(
    title = "D",
    x = "Weeks",
    y = "IRRAT30"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("0" = "black", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("0" = "#b3b3b3", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  geom_point(data = m, aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) #
#+
# theme(legend.position = "none")


# Changing x-axis labels
IRRAT_plot <-
  irrat_plot + scale_x_continuous(
    breaks = c(0, 0.46, 0.99),
    labels = c("0", "24", "52")
  ) +
  scale_y_continuous(limits = c(-1.5, 2)) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.60 (95%CI -1.10 to -0.10), p=0.055",
    x = 0.5, y = 1.8, size = 4
  )


# Get slopes from model

### Get the slopes


# Create list object into which to place model results
slopes <- list()

# Acute slopes
slopes[[1]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes[[2]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes[[3]] <- get_slope(
  .model_obj = fit_IRRAT30, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)

model_irrat <- slopes %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# Make table with flextable


slope_irrat<- qflextable(model_irrat) %>%
 set_table_properties(layout = "autofit") %>%
 set_header_labels(name = "Group", result = "IRRAT30 Slope (95%CI)",
                   "p_value" = "P-Value") %>%
 bold(i = NULL, part = "header")


# Print tables ##
slope_irrat %>% print(preview = "pptx")



## IRITD3 ====
fit_iritd3 <- lmer(IRITD3 ~ Felzartamab * time + (time | Patient), m)

summary(fit_iritd3)

check_model(fit_iritd3)
check_normality(fit_iritd3)

# Get predictions
effect_plot <- ggeffects::ggpredict(fit_iritd3, c("time", "Felzartamab"))

# Plot slopes
IRITD3_plot <-
  ggplot(effect_plot, aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  theme_classic() +
  labs(
    title = "E",
    x = "Weeks",
    y = "IRITD3"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("0" = "black", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("0" = "#b3b3b3", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  geom_point(data = m, aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) #+
# theme(legend.position = "none")


# Changing x-axis labels
IRITD3_plot <-
  IRITD3_plot + scale_x_continuous(
    breaks = c(0, 0.46, 1),
    labels = c("0", "24", "52")
  ) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.20 (95%CI -0.37 to -0.02), p=0.041",
    x = 0.5, y = 0.65, size = 4
  )



# Get slopes from model

### Get the slopes


# Create list object into which to place model results
slopes <- list()

# Acute slopes
slopes[[1]] <- get_slope(
  .model_obj = fit_iritd3, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes[[2]] <- get_slope(
  .model_obj = fit_iritd3, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes[[3]] <- get_slope(
  .model_obj = fit_iritd3, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)

model_IRITD3 <- slopes %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  dplyr::select(name, result, p_value, p_adjusted)

# Make table with flextable


# model_IRITD3<- qflextable(model_IRITD3) %>%
# set_table_properties(layout = "autofit") %>%
#  set_header_labels(name = "Group", result = "IRITD3 Slope (95%CI)",
#                  "p_value" = "P-Value") %>%
# bold(i = NULL, part = "header")


# Print tables ##
model_IRITD3 %>% print(preview = "pptx")



## IRITD5 ====
fit_iritd5 <- lmer(IRITD5 ~ Felzartamab * time + (time | Patient), m)

summary(fit_iritd5)

check_model(fit_iritd5)
check_normality(fit_iritd5)

# Get predictions
effect_plot <- ggeffects::ggpredict(fit_iritd5, c("time", "Felzartamab"))

# Plot slopes
IRITD5_plot <-
  ggplot(effect_plot, aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
  theme_classic() +
  labs(
    title = "F",
    x = "Weeks",
    y = "IRITD5"
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("0" = "black", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("0" = "#b3b3b3", "1" = "#bc3c29"),
    labels = c("Placebo", "Felzartamab")
  ) +
  geom_point(data = m, aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient), alpha = 0, size = 0.1) +
  geom_line(
    data = m, aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient), stat = "smooth", method = "lm",
    alpha = 0.3
  ) #+
# theme(legend.position = "none")


# Changing x-axis labels
IRITD5_plot <- IRITD5_plot + scale_x_continuous(
  breaks = c(0, 0.46, 0.99),
  labels = c("0", "24", "52")
) +
  scale_y_continuous(limits = c(-0.20, 0.8)) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", vjust = -1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  annotate(
    geom = "text",
    label = "Slope difference: -0.23 (95%CI -0.41 to -0.05), p=0.042",
    x = 0.5, y = 0.75, size = 4
  )



# Get slopes from model

### Get the slopes


# Create list object into which to place model results
slopes <- list()

# Acute slopes
slopes[[1]] <- get_slope(
  .model_obj = fit_iritd5, .name = "Placebo",
  # Coefficients: `time`
  .contrasts = c(0, 0, 1, 0) #
)
slopes[[2]] <- get_slope(
  .model_obj = fit_iritd5, .name = "Felzartamab",
  # Coefficients: `time` + `tim:Flzrtmb`
  .contrasts = c(0, 0, 1, 1)
)
slopes[[3]] <- get_slope(
  .model_obj = fit_iritd5, .name = "Intergroup Difference",
  # Coefficients: `tim:Flzrtmb`
  .contrasts = c(0, 0, 0, 1)
)

model_IRITD5 <- slopes %>%
  bind_rows() %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  dplyr::select(name, result, p_value, p_adjusted)


# Make table with flextable


# model_IRITD5<- qflextable(model_IRITD5) %>%
# set_table_properties(layout = "autofit") %>%
# set_header_labels(name = "Group", result = "IRITD5 Slope (95%CI)",
#                   "p_value" = "P-Value") %>%
# bold(i = NULL, part = "header")


# Print tables ##
model_IRITD5 %>% print(preview = "pptx")


# Combine the three slope tables

injury_slopes <- rbind(model_irrat, model_IRITD3, model_IRITD5) %>%
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


# Print tables ##
injury_slopes %>% print(preview = "pptx")


##Save plot and arrange with violin plots =====
# ggarrange(panel_pairs_injury,
          # ggarrange(IRRAT_plot, IRITD3_plot, IRITD5_plot,
                    # nrow = 1, common.legend = T), nrow = 2, ncol = 1)

### Save Tiff

# ggsave(here::here("graphics", "injury_slopes.tiff"),
#        width = 15, height = 9,
#        dpi = 300, compression = "lzw", type="cairo" )

