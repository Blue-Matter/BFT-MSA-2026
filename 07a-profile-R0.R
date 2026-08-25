
library(multiSA)
library(tidyverse)

fit <- readRDS(file.path("model_output", "fit_seasonal_selprior5_08.19d.rds"))

# Single profile of WBFT R0
tictoc::tic()
p <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = seq(100, 350, 25),
  return_models = TRUE,
  cores = 6
)
tictoc::toc()
saveRDS(p, file = file.path("profile", "WBFT_R0_seasonal_selprior5_08.19d.rds"))


# Make figure
if (FALSE) {

  plot_profile <- function(p) {
    MLE <- attr(p$profile, "fitted")[1, 1]

    names_df <- data.frame(
      variable = c("loglike", "loglike_CAL_ymfr", "loglike_I_ymi",
                   "loglike_SC_ymafr",
                   "loglike_tag_mov_ymars", "loglike_CKMR",
                   "logprior_mov",
                   "logprior_sel",
                   "logprior_rdev_ys", "penalty", "objective"),
      variable2 = c("All likelihoods", "Like: Length composition", "Like: Indices",
                    "Like: SOO",
                    "Like: Tags", "Like: CKMR",
                    "Pr: Movement",
                    "Pr: Selectivity",
                    "Pr: Rec devs", "Penalty", "Objective")
    )

    prof_df <- p$profile
    prof_df$loglike_CKMR <- sapply(
      p$fits,
      function(i) {
        ni <- dim(i@report$loglike_I_ymi)[3]
        sum(i@report$loglike_I_ymi[, , ni])
      }
    )
    prof_df$logprior_mov <- sapply(p$fits, function(i) {
      ind <- grepl("mov", names(i@report$logprior_par))
      sum(i@report$logprior_par[ind])
    })
    prof_df$logprior_sel <- sapply(p$fits, function(i) {
      ind <- grepl("sel_pf", names(i@report$logprior_par))
      sum(i@report$logprior_par[ind])
    })

    g <- reshape2::melt(prof_df, id.vars = c("R0_s[2]")) %>%
      mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
      mutate(value = value - min(value), .by = variable) %>%
      left_join(names_df, by = "variable") %>%
      mutate(variable2 = factor(variable2, names_df$variable2)) %>%
      filter(!is.na(variable2)) %>%
      ggplot(aes(`R0_s[2]`, value)) +
      facet_wrap(vars(variable2), scales = "free_y") +
      geom_point() +
      geom_vline(xintercept = MLE, linetype = 2) +
      geom_line() +
      labs(x = expression("WBFT"~~R[0]))
    g
  }

  p <- readRDS(file = file.path("profile", "WBFT_R0_seasonal_selprior5_08.19d.rds"))

  g <- plot_profile(p)
  ggsave("figures/profile/WR0_components.png", g, width = 7, height = 7)
}
