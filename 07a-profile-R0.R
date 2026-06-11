
library(multiSA)
library(tidyverse)

fit <- readRDS(file.path("model_output", "fit_reference1_06.03.2026.rds"))

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
