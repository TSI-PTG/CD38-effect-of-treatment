complex_pivot <- function(data, target, vars, names_to = "Group", values_to = "value") {
    require(tidyverse)
    data <- data %>%
        tidyr::expand_grid(target = target) %>%
        tidyr::nest(.by = target) %>%
        dplyr::mutate(
            data = purrr::pmap(
                list(target, data),
                function(target, data) {
                    target_name <- target %>%
                        stringr::str_remove("_1208Set") %>%
                        stringr::str_remove("^_")
                    data %>%
                        dplyr::select(
                            dplyr::contains(target),
                            dplyr::all_of(vars)
                        ) %>%
                        dplyr::mutate_at(dplyr::vars(dplyr::contains(target)), ~ as.numeric(.) %>% suppressWarnings()) %>%
                        tidyr::pivot_longer(
                            cols = dplyr::contains(target),
                            names_to = names_to,
                            names_pattern = "^(\\w+)_.*",
                            values_to = target_name
                        ) %>%
                        dplyr::mutate_at(dplyr::vars(all_of(names_to)), . %>% stringr::str_extract("[^_]+"))
                }
            )
        )
    purrr::reduce(
        data$data,
        left_join,
        by = c(vars, names_to)
    )
}

