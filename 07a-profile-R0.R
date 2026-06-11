
library(multiSA)
library(tidyverse)

fit <- readRDS(file.path("jitter", "newfit_reference1_06.03.2026.rds"))

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
  plot(p, xlab = "WBFT R0", component = "fn", main = "Change in objective function")

  MLE <- attr(p$profile, "fitted")[1, 1]

  g <- reshape2::melt(p$profile, id.vars = c("R0_s[2]")) %>%
    filter(!variable %in% c("loglike_Cinit_mfr", "fn")) %>%
    mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
    mutate(value = value - min(value), .by = variable) %>%
    ggplot(aes(`R0_s[2]`, value)) +
    facet_wrap(vars(variable), scales = "free_y") +
    geom_point() +
    geom_vline(xintercept = MLE, linetype = 2) +
    geom_line() +
    labs(x = "WBFT R0")

}
