
library(multiSA)
library(tidyverse)

fit <- readRDS(file.path("model_output", "fit_Wprior1_06.22.2026.rds"))

# Single profile of WBFT R0 - SOO1
tictoc::tic()
p <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = c(100, 150, 300, 500, 700, 1000, 1200, 1500, 2000, 3000),
  return_models = TRUE,
  cores = 4
)
tictoc::toc()
saveRDS(p, file = file.path("profile", "Wprior1_WBFT_R0_06.22.2026.rds"))

fit <- readRDS(file.path("model_output", "fit_Wprior2_06.22.2026.rds"))

# Single profile of WBFT R0 - SOO2
tictoc::tic()
p <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = c(100, 150, 300, 500, 700, 1000, 1200, 1500, 2000, 3000),
  return_models = TRUE,
  cores = 4
)
tictoc::toc()
saveRDS(p, file = file.path("profile", "Wprior2_WBFT_R0_06.22.2026.rds"))

# Make figure
if (FALSE) {

  plot_profile <- function(p) {
    MLE <- attr(p$profile, "fitted")[1, 1]

    names_df <- data.frame(
      variable = c("loglike", "loglike_CAL_ymfr", "loglike_I_ymi",
                   #"loglike_SC_ymafr",
                   "loglike_SC_genetic", "loglike_SC_otolith",
                   "loglike_tag_mov_ymars", "logprior",
                   "logprior_spatial", "logprior_sel", "logprior_SSB",
                   #"logprior_par",
                   "logprior_rdev_ys", "penalty", "objective"),
      variable2 = c("All likelihoods", "Like: Length composition", "Like: Indices",
                    #"Like: SOO",
                    "Like: SOO genetic", "Like: SOO otolith",
                    "Like: Tags", "All priors",
                    #"Pr: spatial + sel + SSB prior",
                    "Pr: Spatial", "Pr: Selectivity", "Pr: SSB",
                    "Pr: Rec devs", "Penalty", "Objective")
    )

    prof_df <- p$profile
    prof_df$loglike_SC_otolith <- sapply(
      p$fits,
      function(i) sum(i@report$loglike_SC_ymafr[, , , 1, ])
    )
    prof_df$loglike_SC_genetic <- sapply(
      p$fits,
      function(i) sum(i@report$loglike_SC_ymafr[, , , 2, ])
    )
    prof_df$logprior_spatial <- sapply(p$fits, function(i) {
      ind <- grepl("calc_eqdist", names(i@report$logprior_par))
      sum(i@report$logprior_par[ind])
    })
    prof_df$logprior_sel <- sapply(p$fits, function(i) {
      ind <- grepl("p$sel", names(i@report$logprior_par))
      sum(i@report$logprior_par[ind])
    })
    prof_df$logprior_SSB <- sapply(p$fits, function(i) {
      ind <- grepl("S_yrs", names(i@report$logprior_par))
      sum(i@report$logprior_par[ind])
    })

    g <- reshape2::melt(prof_df, id.vars = c("R0_s[2]")) %>%
      mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
      mutate(value = value - min(value), .by = variable) %>%
      left_join(names_df, by = "variable") %>%
      mutate(variable2 = factor(variable2, names_df$variable2)) %>%
      filter(!is.na(variable2)) %>%
      ggplot(aes(`R0_s[2]`, value)) +
      facet_wrap(vars(variable2), scales = "fixed") +
      geom_point() +
      #coord_transform(x = "log") +
      geom_vline(xintercept = MLE, linetype = 2) +
      geom_line() +
      labs(x = "WBFT R0")
    g
  }

  p <- readRDS(file = file.path("profile", "Wprior1_WBFT_R0_06.22.2026.rds"))

  g <- plot_profile(p) +
    ggtitle("SOO1 + SSB prior")
  ggsave("figures/profile/WR0_SOO1_components.png", g, width = 7, height = 7)

  p <- readRDS(file = file.path("profile", "Wprior2_WBFT_R0_06.22.2026.rds"))
  g <- plot_profile(p) +
    ggtitle("SOO2 + SSB prior")
  ggsave("figures/profile/WR0_SOO2_components.png", g, width = 7, height = 7)

}
