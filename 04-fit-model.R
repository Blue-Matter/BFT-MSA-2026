

library(multiSA)
library(tidyverse)
library(tictoc)
library(parallel)

# Data frame to describe multiple model runs ----
Design <- data.frame(
  initC_scalar = c(0.5, 0.5, 1), # Relative to first year catch
  SSB_prior = c(FALSE, TRUE, FALSE),
  output_name = c("reference_05.20.2026", "Wprior_05.20.2026", "highinitC_05.20.2026"),
  model_name = c("Reference", "WBFT SSB prior", "High eq. catch")
)
readr::write_csv(Design, "tables/Design_05.20.2026.csv")

# Wrapper function that will fit a model for each row in the Design data frame ----
wrapper_fn <- function(x = 1, Design) {

  require(multiSA)
  require(tidyverse)

  dir_save <- "model_input/05.20.2026"

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

  # With WBFT SSB prior
  if (Design$SSB_prior[x]) {
    dat@Dmodel@prior <- c(
      dat@Dmodel@prior,
      paste0("dnorm(log(sum(S_yrs[", match(2018, dat@Dlabel@year), ", , 2])), log(22000), 0.01, log = TRUE)")
    )
  }

  # Rescale equilibrium catch
  dat@Dfishery@Cinit_mfr <- array(
    Design$initC_scalar[x] * dat@Dfishery@Cobs_ymfr[1, , , ],
    c(dat@Dmodel@nm, dat@Dfishery@nf, dat@Dmodel@nr)
  )

  # Check data object
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

  # Movement
  # This can't be automated because there are data for 3 age classes in EBFT but only one age class in WBFT
  # We create only one movement matrix for each stock

  # Estimate EBFT attractivity terms for WATL, EATL, MED
  # Estimate WBFT attractivity terms for GOM, WATL, EATL
  # Attractivity is relative, estimate two parameters for 3 areas and use softmax transformation
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

  # Estimate EBFT and WBFT viscosity term by season (resistance to move from current area)
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

  #### Add selectivity prior?
  add_sel_prior <- TRUE
  if (add_sel_prior) {
    #### Make fishery and rec dist. sel priors ----
    prior_sel <- lapply(1:dat@Dfishery@nf, function(f) {

      # Uninformative prior for length of full selectivity
      p1 <- paste0("dnorm(p$sel_pf[1, ", f, "], 0, 1.5, log = TRUE)")

      #x <- rnorm(1e5, 0, 1.5)
      #hist(plogis(x))

      # Ascending limb with lognormal SD = 0.5
      start_p2 <- round(pars$p$sel_pf[2, f], 2)
      p2 <- paste0("dnorm(p$sel_pf[2, ", f, "], ", start_p2, ", 0.5, log = TRUE)")

      #x <- rnorm(1e5, 0, 0.5)
      #hist(exp(x))

      # Descending limb with lognormal SD = 0.5
      if (grepl("dome", dat@Dfishery@sel_f[f])) {
        start_p3 <- round(pars$p$sel_pf[3, f], 2)
        p3 <- paste0("dnorm(p$sel_pf[3, ", f, "], ", start_p3, ", 0.5, log = TRUE)")
      } else {
        p3 <- NULL
      }

      c(p1, p2, p3)
    }) %>%
      unlist()

    dat@Dmodel@prior <- c(dat@Dmodel@prior, prior_sel)
  }

  #mov <- conv_mov(pars$p$mov_x_marrs[1, , , , 1], pars$p$mov_g_ymars[1, 1, , , 1], pars$p$mov_v_ymas[1, 1, , 1])
  #mov <- conv_mov(x, pars$p$mov_g_ymars[1, 1, , , 2], pars$p$mov_v_ymas[1, 1, , 2])
  #mov[1, , ]

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

  file_out <- paste0("fit_", Design$model_name[x], ".rds")
  saveRDS(fit, file.path("model_output", file_out))

  #fit <- readRDS(file.path("model_output", file_out))
  report(fit, dir = "model_output", filename = paste0("report_", Design$model_name[x]))

  return(invisible(fit))
}


# Fit all models in parallel ----
cl <- parallel::makeCluster(nrow(Design))

tictoc::tic()
fits <- parallel::parLapply(
  cl,
  X = 1:nrow(Design),
  wrapper_fn,
  Design = Design
)
tictoc::toc()

parallel::stopCluster(cl)



