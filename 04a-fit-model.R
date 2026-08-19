

library(multiSA)
library(tidyverse)
library(tictoc)
library(parallel)

source("99-functions-fit.R")

#### Data frame to describe multiple model runs (unique configuration in each row) ----
#### Possible arguments (columns in the data frame)
# input_dir - directory of input model files
# annual - logical for annual or seasonal model (need the correct corresponding directory name)
# movement - logical for estimating movement
# rec_devvector - logical for estimating rec devs as a dev vector (penalty to sum to zero)
# initC_scalar - equilibrium catch (proportion to first year catch)
# lambda_CAL - Lambda factor for length composition
# lambda_SC - Lambda factor for stock mixing
# lambda_tag - Lambda factor for tags (only for seasonal models)
# SC_set - Stock mixing dataset to use (1, 2, or 3)
#          1 = from A. Hanke by year, area, season (April 2026)
#          2 = from I. Fraile by year, area, season (June 2026)
#          3 = from A. Hanke by year and fleet in WATL (genetics only) (August 2026)
# SC_subset - Subset of dataset to fit, only for SC_set = 1 or 2, either "all", "otolith", or "genetic"
# SSB_prior - logical, whether to use CKMR estimate of Western SSB in 2018
# SSB_sd - Numeric, decrease the SE of SSB prior
# spat_prior - logical, whether to use spatial mixing priors
# sel_prior - logical, whether to use sel prior
# Wareas - Either 2 or 3 for number of areas where WBFT can inhabit
# Eareas - Either 2 or 3 for number of areas where EBFT can inhabit
# output_name - Filename to save model
# model_name - Model name for figures

# Annual model
Design <- data.frame(
  input_dir = c("model_input/06.30.2026_annual_2area", "model_input/06.30.2026_annual_4area",
                "model_input/06.30.2026_annual_2area", "model_input/06.30.2026_annual_4area"),
  annual = TRUE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = 0,
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = FALSE,
  SSB_sd = 0.18,
  spat_prior = FALSE,
  sel_prior = c(TRUE, TRUE, FALSE, FALSE),
  fix_sel = c(FALSE, FALSE, TRUE, TRUE),
  Wareas = c(1, 2, 1, 2),
  Eareas = c(1, 2, 1, 2),
  output_name = c("annual_selprior_2area_08.19", "annual_selprior_4area_08.19",
                  "annual_fixsel_2area_08.19", "annual_fixsel_4area_08.19"),
  model_name = c("Annual, 2 area", "Annual, 4 area",
                 "Annual, 2 area, fixsel", "Annual, 4 area, fixsel")
)
readr::write_csv(Design, "tables/Design_08.19.2026_annual.csv")

# Fit all models in parallel or in a loop ----
do_parallel <- TRUE

if (do_parallel) {

  cl <- parallel::makeCluster(4) # Need to limit the number of cores due to memory constraints
  tictoc::tic()
  fits <- parallel::parLapplyLB(
    cl,
    X = 1:nrow(Design),
    wrapper_fn,
    Design = Design
  )
  tictoc::toc()

  parallel::stopCluster(cl)

} else {

  fits <- list()

  for (i in 1:nrow(Design)) {
    message("Model ", i)
    tictoc::tic()
    fits[[i]] <- wrapper_fn(i, Design = Design)
    tictoc::toc()
  }

}

