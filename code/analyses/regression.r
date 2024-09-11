# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(lme4) # install.packages("lme4")
library(lmerTest) # install.packages("lmerTest")
library(performance) # install.packages("performance")
library(ggeffects) # install.packages("ggeffects")
library(flextable) # install.packages("flextable")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("code/functions/get_slope_function.r")
# load reference set
load("data/cd38.RData")



# WRANGLE THE DATA ####
data <- set %>%
    pData() %>%
    dplyr::select(Center, Patient, Felzartamab, Group, IRRAT30, IRITD3, IRITD5) %>%
    dplyr::filter(Patient %nin% c(15, 18)) %>%
    mutate(
        time = case_when(
            Group == "Index" ~ 0,
            Group == "FU1" ~ 0.46,
            Group == "FU2" ~ 0.99
        ),
        linetype = "solid" %>% factor()
    )


# FIT MODELS ####
# IRRAT30
fit_IRRAT30 <- lmer(IRRAT30 ~ Felzartamab * time + (time | Patient), data)
# test model assumptions
performance::check_model(fit_IRRAT30)
performance::check_normality(fit_IRRAT30)
# test linearity
fit_IRRAT30_nonlinear <- lmer(IRRAT30 ~ Felzartamab * poly(time, 2) + (time | Patient), data)
anova(fit_IRRAT30, fit_IRRAT30_nonlinear)
# Get predictions
effect_plot_IRRAT30 <- ggeffects::ggpredict(fit_IRRAT30, c("time", "Felzartamab"))

# IRITD3
fit_IRITD3 <- lmer(IRITD3 ~ Felzartamab * time + (time | Patient), data)
# test model assumptions
performance::check_model(fit_IRITD3)
performance::check_normality(fit_IRITD3)
# test linearity
fit_IRITD3_nonlinear <- lmer(IRITD3 ~ Felzartamab * poly(time, 2) + (time | Patient), data)
anova(fit_IRITD3, fit_IRITD3_nonlinear)
# Get predictions
effect_plot_IRITD3 <- ggeffects::ggpredict(fit_IRITD3, c("time", "Felzartamab"))

# IRITD5
fit_IRITD5 <- lmer(IRITD5 ~ Felzartamab * time + (time | Patient), data)
# test model assumptions
performance::check_model(fit_IRITD5)
performance::check_normality(fit_IRITD5)
# test linearity
fit_IRITD5_nonlinear <- lmer(IRITD5 ~ Felzartamab * poly(time, 2) + (time | Patient), data)
anova(fit_IRITD5, fit_IRITD5_nonlinear)
# Get predictions
effect_plot_IRITD5 <- ggeffects::ggpredict(fit_IRITD5, c("time", "Felzartamab"))


# EXTRACT THE SLOPES FROM THE MODEL FITS ####
# IRRAT30
# Create list object into which to place model results
slopes_IRRAT30 <- list()
# Acute slopes
slopes_IRRAT30[[1]] <- get_slope(
    .model_obj = fit_IRRAT30, .name = "Placebo",
    .contrasts = c(0, 0, 1, 0) #
)
slopes_IRRAT30[[2]] <- get_slope(
    .model_obj = fit_IRRAT30, .name = "Felzartamab",
    .contrasts = c(0, 0, 1, 1)
)
slopes_IRRAT30[[3]] <- get_slope(
    .model_obj = fit_IRRAT30, .name = "Intergroup Difference",
    .contrasts = c(0, 0, 0, 1)
)
model_IRRAT30 <- slopes_IRRAT30 %>%
    bind_rows() %>%
    mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
    dplyr::select(name, result, p_value, p_adjusted)

# IRITD3
# Create list object into which to place model results
slopes_IRITD3 <- list()
# Acute slopes
slopes_IRITD3[[1]] <- get_slope(
    .model_obj = fit_IRITD3, .name = "Placebo",
    .contrasts = c(0, 0, 1, 0) #
)
slopes_IRITD3[[2]] <- get_slope(
    .model_obj = fit_IRITD3, .name = "Felzartamab",
    .contrasts = c(0, 0, 1, 1)
)
slopes_IRITD3[[3]] <- get_slope(
    .model_obj = fit_IRITD3, .name = "Intergroup Difference",
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
    .contrasts = c(0, 0, 1, 0) #
)
slopes_IRITD5[[2]] <- get_slope(
    .model_obj = fit_IRITD5, .name = "Felzartamab",
    .contrasts = c(0, 0, 1, 1)
)
slopes_IRITD5[[3]] <- get_slope(
    .model_obj = fit_IRITD5, .name = "Intergroup Difference",
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


# WRANGLE RESULTS FOR PLOTTING ####
regression_injury <- tibble(
    variable = c("IRRAT30", "IRITD3", "IRITD5"),
    reference_data = data %>% list(),
    plot_data = list(effect_plot_IRRAT30, effect_plot_IRITD3, effect_plot_IRITD5),
    slope_tables = list(table_slope_IRRAT30, table_slope_IRITD3, table_slope_IRITD5)
)


# EXPORT RESULTS FOR PLOTTING ####
saveDir <- "results/"
save(regression_injury, file = paste(saveDir, "regression_injury.RData", sep = ""))
