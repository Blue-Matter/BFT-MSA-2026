

library(multiSA)
library(tidyverse)

Design <- readr::read_csv("tables/Design_08.19.2026_annual.csv")
model_name <- paste0(Design$output_name[1:4], ".rds")

for (j in 1:length(model_name)) {
  i <- model_name[j]
  fit <- readRDS(file.path("model_output", paste0("fit_", i)))

  tictoc::tic()
  fits <- do_jitter(fit, n = 10, amount = 0.25, seed = 23, cores = 5, use_fitted = FALSE, return_models = TRUE, do_sd = FALSE)
  tictoc::toc()

  saveRDS(fits, file = file.path("jitter", paste0("jitter_", i)))

  # Need time to compute max gradient
  like <- do.call(rbind, lapply(fits, get_likelihood_components))
  saveRDS(like, file = file.path("jitter", paste0("jitter_loglik_", i)))
}

if (FALSE) {
  j <- 1
  i <- Design$output_name[j]
  like <- readRDS(file = file.path("jitter", paste0("jitter_loglik_", i, ".rds")))

  names_df <- data.frame(
    variable = c("loglike", "loglike_CAL_ymfr", "loglike_I_ymi", "loglike_SC_ymafr",
                 "loglike_tag_mov_ymars", "logprior", "logprior_par",
                 "logprior_rdev_ys", "penalty", "objective", "maxgrad"),
    variable2 = c("All likelihoods", "Like: Length composition", "Like: Indices", "Like: SOO",
                  "Like: Tags", "All priors", "Pr: spatial prior",
                  "Pr: Rec devs", "Penalty", "Objective", "Max. Gradient")
  )

  g <- like %>%
    mutate(Jitter = 1:n()) %>%
    reshape2::melt(id.vars = "Jitter") %>%
    filter(!variable %in% c("loglike_Cinit_mfr", "conv", "fn")) %>%
    mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
    #mutate(value = value - min(value), .by = variable) %>%
    left_join(names_df) %>%
    mutate(variable2 = factor(variable2, names_df$variable2)) %>%
    ggplot(aes(Jitter, value)) +
    facet_wrap(vars(variable2), scales = "free_y") +
    geom_point() +
    geom_line() +
    labs(x = "Jitter Run", title = Design$model_name[j])
  ggsave("figures/jitter_run.png", g, height = 6, width = 8)

}
