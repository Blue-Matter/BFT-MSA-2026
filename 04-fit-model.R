

library(multiSA)
library(tidyverse)

dir_save <- "model_input/04.30.2026"

# Make MSA data object from saved objects
dat <- new(
  "MSAdata",
  Dmodel = readRDS(file.path(dir_save, "Dmodel.rds")),
  Dstock = readRDS(file.path(dir_save, "Dstock_A.rds")),
  Dfishery = readRDS(file.path(dir_save, "Dfishery_SOO_genetic.rds")),
  Dsurvey = readRDS(file.path(dir_save, "Dsurvey.rds")),
  Dtag = readRDS(file.path(dir_save, "Dtag.rds")),
  Dlabel = readRDS(file.path(dir_save, "Dlabel.rds"))
)

dat <- check_data(dat)

# EBFT recruits in MED
# WBFT recruits in GOM, WATL
log_recdist_rs <- matrix(0, dat@Dmodel@nr, dat@Dmodel@ns)
log_recdist_rs[-4, 1] <- -1000
log_recdist_rs[3:4, 2] <- -1000

# We only estimate one recruitment distribution parameter: proportion recruitment in GOM for WBFT
# (Fix WATL parameter below)
map_recdist_rs <- matrix(NA, dat@Dmodel@nr, dat@Dmodel@ns)
map_recdist_rs[2, 2] <- 1


parameters_start <- list(
  log_recdist_rs = log_recdist_rs,
  R0_s = c(20000, 5000),
  h_s = c(0.99, 0.6),
  log_sdr_s = rep(log(0.5), 2),
  log_initF_mfr = array(log(0.01), c(dat@Dmodel@nm, dat@Dfishery@nf, dat@Dmodel@nr))
)
pars <- make_parameters(
  dat,
  start = parameters_start,
  map = list(log_recdist_rs = factor(map_recdist_rs)),
  est_mov = "gravity_fixed"
)


tictoc::tic()
fit <- fit_MSA(
  dat,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()
saveRDS(fit, file = "model_output/fit_04.30.2026.rds")

fit <- readRDS("model_output/fit_04.30.2026.rds")
report(fit, dir = "model_output", filename = "report_04.30.2026")
