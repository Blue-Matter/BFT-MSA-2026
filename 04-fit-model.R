

library(multiSA)
library(tidyverse)

dir_save <- "model_input/04.30.2026"

#### Make MSA data object from saved objects ----
dat <- new(
  "MSAdata",
  Dmodel = readRDS(file.path(dir_save, "Dmodel.rds")),
  Dstock = readRDS(file.path(dir_save, "Dstock_A.rds")),
  Dfishery = readRDS(file.path(dir_save, "Dfishery.rds")),
  Dsurvey = readRDS(file.path(dir_save, "Dsurvey.rds")),
  Dtag = readRDS(file.path(dir_save, "Dtag.rds")),
  Dlabel = readRDS(file.path(dir_save, "Dlabel.rds"))
)
#dat@Dmodel@prior <- dat@Dmodel@prior[-c(1:2)]
dat@Dfishery@fcomp_like <- "multinomial"

dat <- check_data(dat)

#### Starting parameters ----

# EBFT recruits in MED
# WBFT recruits in GOM, WATL
log_recdist_rs <- matrix(0, dat@Dmodel@nr, dat@Dmodel@ns)
log_recdist_rs[-4, 1] <- -1000
log_recdist_rs[3:4, 2] <- -1000

parameters_start <- list(
  log_recdist_rs = log_recdist_rs,
  R0_s = c(20000, 10000),
  h_s = c(0.99, 0.6),
  log_sdr_s = log(c(0.5, 0.5))
)

#### Fixing parameters ----
# NA = fix, integer = estimate, shared integers = shared parameter value

# We only estimate one recruitment distribution parameter: proportion recruitment in GOM for WBFT
# (Fix WATL parameter below)
map_recdist_rs <- matrix(NA, dat@Dmodel@nr, dat@Dmodel@ns)
map_recdist_rs[1, 2] <- 1

# Recruitment deviations
# All years for EBFT
# Only after 1960 for WBFT
map_log_rdev_ys <- matrix(NA, dat@Dmodel@ny, dat@Dmodel@ns)
map_log_rdev_ys[, 1] <- map_log_rdev_ys[dat@Dlabel@year > 1960, 2] <- TRUE
map_log_rdev_ys[!is.na(map_log_rdev_ys)] <- 1:sum(map_log_rdev_ys, na.rm = TRUE)


map_g_ymars <- array(c(dat@Dmodel@ny, dat@Dmodel@nm, dat@Dmodel@na, dat@Dmodel@nr, dat@Dmodel@ns))
# Movement
#
# This can't be automated because there are data for 3 age classes in EBFT but only one age class in WBFT
# We create only one movement matrix for each stock

# Estimate EBFT gravity term (to stay in current area) for all areas except GOM
# Estimate WBFT gravity term (to stay in current area) for all areas except MED
# Softmax transformation: remove a second area
map_g_ymars <- array(NA, c(dat@Dmodel@ny, dat@Dmodel@nm, dat@Dmodel@na, dat@Dmodel@nr, dat@Dmodel@ns))
for (s in 1:dat@Dmodel@ns) {
  g_s <- matrix(NA_real_, dat@Dmodel@nm, dat@Dmodel@nr)
  if (s == 1) {
    g_s[, -c(1:2)] <- TRUE
    g_s[, -c(1:2)] <- 1:sum(g_s, na.rm = TRUE)
  } else {
    g_s[, -c(3:4)] <- TRUE
    g_s[, -c(3:4)] <- max(map_g_ymars, na.rm = TRUE) + 1:sum(g_s, na.rm = TRUE)
  }
  map_g_ymars[, , , , s] <- array(g_s, c(dat@Dmodel@nm, dat@Dmodel@nr, dat@Dmodel@ny, dat@Dmodel@na)) %>%
    aperm(c(3, 1, 4, 2))
}
range(map_g_ymars, na.rm = TRUE)

# Estimate EBFT viscosity term (resistance to move) for all areas except GOM
# Estimate WBFT viscosity term (resistance to move) for all areas except MED
# Softmax transformation: remove a second area
map_v_ymas <- matrix(seq(1, dat@Dmodel@nm * dat@Dmodel@ns), dat@Dmodel@nm, dat@Dmodel@ns) %>%
  array(c(dat@Dmodel@nm, dat@Dmodel@ns, dat@Dmodel@ny, dat@Dmodel@na)) %>%
  aperm(c(3, 1, 4, 2))

map <- list(
  log_recdist_rs = factor(map_recdist_rs),
  log_rdev_ys = factor(map_log_rdev_ys),
  mov_g_ymars = map_g_ymars,
  mov_v_ymas = map_v_ymas
)

#### Make full parameter and map listss ----
pars <- make_parameters(
  dat,
  start = parameters_start,
  map = map,
  est_mov = "gravity_fixed"
)

# Manually check movement matrix setup
#x <- array(0, c(36, 4, 4))
#x[, , 1] <- -1000
#x[, 1, ] <- 1000

#mov <- conv_mov(pars$p$mov_x_marrs[1, , , , 1], pars$p$mov_g_ymars[1, 1, , , 1], pars$p$mov_v_ymas[1, 1, , 1])
#mov <- conv_mov(x, pars$p$mov_g_ymars[1, 1, , , 2], pars$p$mov_v_ymas[1, 1, , 2])
#mov[1, , ]

tictoc::tic()
fit <- fit_MSA(
  dat,
  pars$p,
  pars$map,
  pars$random,
  run_model = FALSE,
  do_sd = TRUE
)
tictoc::toc()
saveRDS(fit, file = "model_output/fit_04.30.2026.rds")

fit <- readRDS("model_output/fit_04.30.2026.rds")
report(fit, dir = "model_output", filename = "report_04.30.2026")
