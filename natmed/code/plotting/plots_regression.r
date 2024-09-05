# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(ggpubr) # install.packages("ggpubr")
# load regression results
loadDir <- "natmed/results/"
load(paste(loadDir, "regression_injury.RData", sep = ""))


# MAKE PLOTS ####
plots_regression <- regression_injury %>%
    mutate(
        plot = pmap(
            list(variable, reference_data, plot_data),
            function(variable, reference_data, plot_data) {
                label <- case_when(
                    variable == "IRRAT30" ~ "Injury-repair associated\n(IRRAT30)",
                    variable == "IRITD3" ~ "Injury-repair induced, day 3\n(IRITD3)",
                    variable == "IRITD5" ~ "Injury-repair induced, day 5\n(IRITD5)"
                )
                ylim <- case_when(
                    variable == "IRRAT30" ~ c(-1.5, 2),
                    variable == "IRITD3" ~ c(-0.7, 0.7),
                    variable == "IRITD5" ~ c(-0.20, 0.8)
                )
                plot_data %>%
                    ggplot(aes(x = x, y = predicted, color = group)) +
                    geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
                    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
                    geom_point(
                        data = reference_data,
                        aes_string(x = "time", y = variable, color = factor("Felzartamab"), group = "Patient"),
                        alpha = 0,
                        size = 0.1
                    ) +
                    geom_line(
                        data = reference_data,
                        aes_string(x = "time", y = variable, color = factor("Felzartamab"), group = "Patient", linetype = "linetype"),
                        stat = "smooth",
                        method = "lm",
                        alpha = 0.3,
                        show.legend = FALSE
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
                    scale_y_continuous(limits = ylim) +
                    labs(
                        x = "Time post-treatment (weeks)",
                        y = label
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
            }
        )
    )


# GENERATE GROBS FOR SLOPE TABLES ####
plots_regression <- plots_regression %>%
    mutate(
        grob_table = map(
            slope_tables,
            function(slope_tables) {
                slope_tables %>%
                    gen_grob(fit = "fixed", just = "centre") %>%
                    as_ggplot()
            }
        )
    )


# SAVE THE PLOT DATA ####
saveDir <- "natmed/results/"
save(plots_regression, file = paste(saveDir, "plots_regression.RData", sep = ""))

