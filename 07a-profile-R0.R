
library(multiSA)
library(tidyverse)

fit <- readRDS(file.path("model_output", "newfit_reference1_06.03.2026.rds"))

# Single profile of WBFT R0
tictoc::tic()
p <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = c(100, 150, 300, 500, 700, 1000, 1200, 1500, 2000, 3000),
  return_models = FALSE,
  cores = 4
)
tictoc::toc()
saveRDS(p, file = file.path("profile", "reference1_WBFT_R0_06.03.2026.rds"))

if (FALSE) {
  p <- readRDS(file = file.path("profile", "reference1_WBFT_R0_06.03.2026.rds"))

  png("figures/profile/WR0_ref1.png", height = 4, width = 5, units = "in", res = 400)
  par(mar = c(5, 4, 1, 1))
  plot(p, xlab = "WBFT R0", component = "fn", ylab = "Change in objective function", main = "SOO1")
  dev.off()

  MLE <- attr(p$profile, "fitted")[1, 1]

  names_df <- data.frame(
    variable = c("loglike", "loglike_CAL_ymfr", "loglike_I_ymi", "loglike_SC_ymafr",
                 "loglike_tag_mov_ymars", "logprior", "logprior_par",
                 "logprior_rdev_ys", "penalty", "objective"),
    variable2 = c("All likelihoods", "Like: Length composition", "Like: Indices", "Like: SOO",
                  "Like: Tags", "All priors", "Pr: spatial prior",
                  "Pr: Rec devs", "Penalty", "Objective")
  )

  g <- reshape2::melt(p$profile, id.vars = c("R0_s[2]")) %>%
    filter(!variable %in% c("loglike_Cinit_mfr", "fn")) %>%
    mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
    mutate(value = value - min(value), .by = variable) %>%
    left_join(names_df) %>%
    mutate(variable2 = factor(variable2, names_df$variable2)) %>%
    ggplot(aes(`R0_s[2]`, value)) +
    facet_wrap(vars(variable2), scales = "free_y") +
    geom_point() +
    #coord_transform(x = "log") +
    geom_vline(xintercept = MLE, linetype = 2) +
    geom_line() +
    ggtitle("SOO1") +
    labs(x = "WBFT R0")
  ggsave("figures/profile/WR0_ref1_components.png", g, width = 8, height = 6)


  p <- readRDS(file = file.path("profile", "reference2_WBFT_R0_06.03.2026.rds"))

  png("figures/profile/WR0_ref2.png", height = 4, width = 5, units = "in", res = 400)
  par(mar = c(5, 4, 1, 1))
  plot(p, xlab = "WBFT R0", component = "fn", ylab = "Change in objective function", main = "SOO2")
  dev.off()

  MLE <- attr(p$profile, "fitted")[1, 1]
  g <- reshape2::melt(p$profile, id.vars = c("R0_s[2]")) %>%
    filter(!variable %in% c("loglike_Cinit_mfr", "fn")) %>%
    mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
    mutate(value = value - min(value), .by = variable) %>%
    left_join(names_df) %>%
    mutate(variable2 = factor(variable2, names_df$variable2)) %>%
    ggplot(aes(`R0_s[2]`, value)) +
    facet_wrap(vars(variable2), scales = "free_y") +
    geom_point() +
    #coord_transform(x = "log") +
    geom_vline(xintercept = MLE, linetype = 2) +
    geom_line() +
    ggtitle("SOO2") +
    labs(x = "WBFT R0")
  ggsave("figures/profile/WR0_ref2_components.png", g, width = 8, height = 6)

}
