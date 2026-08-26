
# Wrapper function that will fit a model for each row in the Design data frame ----
wrapper_fn <- function(x = 1, Design) {

  require(multiSA)
  require(tidyverse)

  input_dir <- Design$input_dir[x]

  #### Make MSA data object from saved objects ----
  dat <- new(
    "MSAdata",
    Dmodel = readRDS(file.path(input_dir, "Dmodel.rds")),
    Dstock = readRDS(file.path(input_dir, "Dstock_A.rds")),
    Dfishery = readRDS(file.path(input_dir, paste0("Dfishery_SOO", Design$SC_set[x], ".rds"))),
    Dsurvey = readRDS(file.path(input_dir, "Dsurvey.rds")),
    Dtag = readRDS(file.path(input_dir, "Dtag.rds")),
    Dlabel = readRDS(file.path(input_dir, "Dlabel.rds"))
  )

  # Downweight CAL
  dat@Dfishery@lambdaCAL_f <- Design$lambda_CAL[x]

  # Francis weights
  #if (Design$annual[x]) {
  #  wt <- c(0.09, 0.23, 0.33, 0.51, 0.98, 0.19, 0.44, 2.42, 7.10, 0.09, 0.97, 0.72, 0.32, 0.37, 2.06, 5.94, 0.04, 0.10)
  #  for (f in 1:length(wt)) {
  #    if (wt[f] < 1) dat@Dfishery@CALN_ymfr[, , f, ] <- wt[f] * dat@Dfishery@CALN_ymfr[, , f, ]
  #  }
  #}

  # Downweight SC
  dat@Dfishery@lambdaSC_f <- Design$lambda_SC[x]
  #dat@Dfishery@SCstdev_ymafrs[] <- dat@Dfishery@SCstdev_ymafrs + 0.1

  #if (Design$SC_subset[x] == "otolith") dat@Dfishery@lambdaSC_f[2] <- 0 # Set genetic to zero
  #if (Design$SC_subset[x] == "genetic") dat@Dfishery@lambdaSC_f[1] <- 0 # Set otolith to zero

  # Downweight tag
  if (!Design$annual[x]) dat@Dtag@lambdaTag_s <- Design$lambda_tag[x]

  # Rescale equilibrium catch
  dat@Dfishery@Cinit_mfr <- array(
    Design$initC_scalar[x] * dat@Dfishery@Cobs_ymfr[1, , , ],
    c(dat@Dmodel@nm, dat@Dfishery@nf, dat@Dmodel@nr)
  )

  # Check data object
  dat <- check_data(dat)

  # Constrain rec devs as devvector (sums to zero)
  if (Design$rec_devvector[x]) {
    prior_recdev <- paste0("dnorm(sum(p$log_rdev_ys[, ", 1:2, "]), 0, 0.01, log = TRUE)")
    dat@Dmodel@prior <- c(dat@Dmodel@prior, prior_recdev)
  }

  # Add (or remove) spatial prior
  if (!Design$movement[x] || !Design$spat_prior[x]) {
    dat@Dmodel@prior <- dat@Dmodel@prior[!grepl("calc_eqdist", dat@Dmodel@prior)]
  }

  #### Starting parameters ----
  log_recdist_rs <- matrix(0, dat@Dmodel@nr, dat@Dmodel@ns)

  if (Design$movement[x]) {
    # EBFT recruits in MED
    log_recdist_rs[1:3, 1] <- -1000
  } else if (Design$Eareas[x] == 1) {
    log_recdist_rs[1, 1] <- -1000
  } else if (Design$Eareas[x] == 2) {
    # EBFT recruits in EATL/MED
    log_recdist_rs[1:2, 1] <- -1000
  } else if (Design$Eareas[x] == 3) {
    # EBFT recruits in WATL/EATL/MED
    log_recdist_rs[1, 1] <- -1000
  }

  if (Design$Wareas[x] == 1) {
    log_recdist_rs[2, 2] <- -1000
  } else {
    # WBFT recruits in GOM, WATL
    log_recdist_rs[3:4, 2] <- -1000
  }

  parameters_start <- list(
    log_recdist_rs = log_recdist_rs,
    #R0_s = c(5000, 1000),
    R0_s = c(5000, 5000),
    h_s = c(0.99, 0.6),
    log_sdr_s = log(c(0.5, 0.2))
  )

  #### Fixing parameters ----
  # NA = fix, integer = estimate, shared integers = shared parameter value

  # Recruitment distribution parameter: proportion recruitment among areas by stock
  # There is always one less degree of freedom than the number of areas for distributing recruitment
  map_recdist_rs <- matrix(NA, dat@Dmodel@nr, dat@Dmodel@ns)

  if (Design$Wareas[x] == 1) {
    map_recdist_rs[, 2] <- NA
  } else {
    # WBFT always recruit to GOM and WATL
    map_recdist_rs[1, 2] <- 1
  }

  if (Design$movement[x]) {
    # EBFT recruits in MED (don't estimate recruitment distribution)
  } else if (Design$Eareas[x] == 1) {
    map_recdist_rs[, 1] <- NA
  } else if (Design$Eareas[x] == 2) {
    # EBFT recruits in EATL/MED (est proportion in MED)
    map_recdist_rs[4, 1] <- 2
  } else if (Design$Eareas[x] == 3) {
    # EBFT recruits in WATL/EATL/MED (est proportions in EATL/MED)
    map_recdist_rs[3:4, 1] <- 2:3
  }

  # Recruitment deviations
  # Don't estimate last two years
  # Otherwise: All years for EBFT, Only after 1960 for WBFT
  map_log_rdev_ys <- matrix(0, dat@Dmodel@ny, dat@Dmodel@ns)
  map_log_rdev_ys[, 1] <- map_log_rdev_ys[dat@Dlabel@year > 1960, 2] <- 1
  map_log_rdev_ys[seq(-1, 0) + dat@Dmodel@ny, ] <- 0

  #map_log_rdev_ys[, 2] <- 0

  map_log_rdev_ys[map_log_rdev_ys > 0] <- 1:sum(map_log_rdev_ys, na.rm = TRUE)
  map_log_rdev_ys[map_log_rdev_ys == 0] <- NA

  # Stock selectivity
  map_log_q_fs <- matrix(NA, dat@Dfishery@nf, dat@Dmodel@ns)
  if (Design$est_stocksel[x]) {
    map_log_q_fs[colSums(dat@Dfishery@SC_ff) > 0, 2] <- 1
    map_log_q_fs[!is.na(map_log_q_fs)] <- 1:sum(map_log_q_fs, na.rm = TRUE)
  } else {
    map_log_q_fs <- matrix(NA, dat@Dfishery@nf, dat@Dmodel@ns)
  }

  map <- list(
    log_recdist_rs = factor(map_recdist_rs),
    log_rdev_ys = factor(map_log_rdev_ys),
    log_q_fs = factor(map_log_q_fs)
  )

  # Movement
  # This can't be automated because there are data for 3 age classes in EBFT but only one age class in WBFT
  # We create only one movement matrix for each stock

  # Estimate EBFT attractivity terms for WATL, EATL, MED
  # Estimate WBFT attractivity terms for GOM, WATL, EATL
  # Attractivity is relative, there is always one less degree of freedom than the number of areas
  # For example, estimate two parameters for 3 areas and use softmax transformation
  if (Design$movement[x] && !Design$annual[x]) {
    prior_mov <- NULL
    map_g_ymars <- array(NA, c(dat@Dmodel@ny, dat@Dmodel@nm, dat@Dmodel@na, dat@Dmodel@nr, dat@Dmodel@ns))

    for (s in 1:dat@Dmodel@ns) {
      g_s <- matrix(NA_real_, dat@Dmodel@nm, dat@Dmodel@nr)
      if (s == 1) {
        i <- 3:4
        g_s[, i] <- TRUE
        g_s[, i] <- 1:sum(g_s, na.rm = TRUE)
      } else {
        if (Design$Wareas[x] == 2) {
          dat@Dstock@presence_rs[3, 2] <- FALSE # Turn off WATL presence in EATL (r = 3, s = 2)
          i <- 1
        } else if (Design$Wareas[x] == 3) {
          i <- 1:2
        }
        g_s[, i] <- TRUE
        g_s[, i] <- max(map_g_ymars, na.rm = TRUE) + 1:sum(g_s, na.rm = TRUE)
      }
      map_g_ymars[, , , , s] <- array(g_s, c(dat@Dmodel@nm, dat@Dmodel@nr, dat@Dmodel@ny, dat@Dmodel@na)) %>%
        aperm(c(3, 1, 4, 2))

      prior_mov <- c(
        prior_mov,
        sapply(1:dat@Dmodel@nm, function(m) {
          sapply(i, function(r) {
            paste0("dnorm(p$mov_g_ymars[1, ", m, ", 1, ", r, ", ", s, "], 0, 1.5, log = TRUE)")
          })
        }) %>% as.character()
      )
    }
    range(map_g_ymars, na.rm = TRUE)

    # Estimate EBFT and WBFT viscosity term by season (resistance to move from current area)
    map_v_ymas <- matrix(seq(1, dat@Dmodel@nm * dat@Dmodel@ns), dat@Dmodel@nm, dat@Dmodel@ns) %>%
      array(c(dat@Dmodel@nm, dat@Dmodel@ns, dat@Dmodel@ny, dat@Dmodel@na)) %>%
      aperm(c(3, 1, 4, 2))

    prior_mov <- c(
      prior_mov,
      sapply(1:dat@Dmodel@nm, function(m) {
        sapply(1:dat@Dmodel@ns, function(s) {
          paste0("dnorm(p$mov_v_ymas[1, ", m, ", 1, ", s, "], 0, 1.5, log = TRUE)")
        })
      }) %>% as.character()
    )

    map$mov_g_ymars <- map_g_ymars
    map$mov_v_ymas <- map_v_ymas

    # Add logit prior to avoid going to extreme values
    dat@Dmodel@prior <- c(dat@Dmodel@prior, prior_mov)
  }

  #if (Design$fix_sel[x]) {
  #  map$sel_pf <- matrix(c(NA, TRUE, NA), 3, max(dat@Dfishery@sel_block_yf))
  #  for (f in 1:length(dat@Dfishery@sel_f)) {
  #    if (grepl("dome", dat@Dfishery@sel_f[f])) map$sel_pf[3, f] <- TRUE
  #  }
  #  map$sel_pf[!is.na(map$sel_pf)] <- 1:sum(map$sel_pf, na.rm = TRUE)
  #}

  #### Make full parameter and map lists ----
  pars <- make_parameters(
    dat,
    start = parameters_start,
    map = map,
    est_mov = ifelse(Design$movement[x], "gravity_fixed", "none"),
    silent = TRUE
  )

  # Manually check movement matrix setup
  #x <- array(0, c(36, 4, 4))
  #x[, , 1] <- -1000
  #x[, 1, ] <- 1000

  #### Make fishery sel priors ----
  if (Design$sel_prior[x]) {
    nb <- length(dat@Dfishery@sel_f)
    prior_sel <- lapply(1:nb, function(f) {

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

      p3 <- NULL

      c(p1, p2, p3)
    }) %>%
      unlist()

    dat@Dmodel@prior <- c(dat@Dmodel@prior, prior_sel)
  }

  dat@Dsurvey@lambdaI_i <- rep(1, dat@Dsurvey@ni)

  # Add CKMR estimate of WBFT SSB
  if (!Design$SSB_prior[x]) {
    # Set CKMR survey likelihood weight to zero
    dat@Dsurvey@lambdaI_i[dat@Dlabel@index == "WBFT_CKMR"] <- 0
  } else {
    dat@Dsurvey@Isd_ymi[match(2018, dat@Dlabel@year), ifelse(Design$annual[x], 1, 2), dat@Dsurvey@ni] <- Design$SSB_sd[x]
  }

  # Fit model
  fit <- fit_MSA(
    dat,
    pars$p,
    pars$map,
    pars$random,
    run_model = TRUE,
    do_sd = !Design$fix_sel[x]
  )

  if (Design$fix_sel[x]) {

    pars$p$sel_pf[] <- fit@obj$env$parList(fit@opt$par)$sel_pf
    pars$map$sel_pf <- factor(array(NA, dim(pars$p$sel_pf)))

    fit <- fit_MSA(
      dat,
      pars$p,
      pars$map,
      pars$random,
      run_model = TRUE,
      do_sd = TRUE
    )
  }

  # Save model object
  file_out <- paste0("fit_", Design$output_name[x], ".rds")
  saveRDS(fit, file.path("model_output", file_out))

  # Can take a while to render document
  if (FALSE) {
    report(
      fit,
      name = Design$model_name[x],
      dir = "model_output",
      filename = paste0("report_", Design$output_name[x]),
      open_file = FALSE
    )
  }

  return(invisible(fit))
}

