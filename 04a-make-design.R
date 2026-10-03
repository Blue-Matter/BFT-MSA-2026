

library(multiSA)
library(tidyverse)
library(tictoc)
library(parallel)

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
# lambda_SSB - Lambda factor for CKMR estimate of Western SSB in 2018
# spat_prior - logical, whether to use spatial mixing priors
# sel_prior - logical, whether to use sel prior
# Wareas - Either 2 or 3 for number of areas where WBFT can inhabit
# Eareas - Either 2 or 3 for number of areas where EBFT can inhabit
# output_name - Filename to save model
# model_name - Model name for figures

# Seasonal, no movement
# 4 models:
# (1) no CKMR
# (2) With CKMR estimate, SD = 0.18
# (2a) With CKMR estimate, upweight lambda = 10
# (3) 2a, with CKMR + stock mixing in West Atlantic
Design <- data.frame(
  input_dir = "model_input/10.02.2026",
  annual = FALSE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = 0,
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  lambda_SSB = c(0, 1, 10),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = FALSE,
  minI_SD = 0.1,
  minSC_SD = NA,
  Wareas = 2,
  Eareas = 2,
  output_name = paste0("seasonal_selprior", 1:3, "_10.02"),
  model_name = c("(1) NM: no CKMR", "(2) NM: CKMR", "(2a) NM: CKMR, upweight")
)
readr::write_csv(Design, "tables/Design_10.02.2026_seasonal.csv")

Design <- data.frame(
  input_dir = "model_input/06.30.2026",
  annual = FALSE,
  movement = c(FALSE, TRUE),
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = 1,
  lambda_tag = c(0, 1),
  SC_set = 3,
  SC_subset = "all",
  lambda_SSB = 10,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = TRUE,
  minI_SD = 0.1,
  minSC_SD = NA,
  Wareas = 2,
  Eareas = 3,
  output_name = paste0("seasonal_selprior", 4:5, "_10.02"),
  model_name = c("(3) NM: CKMR+mixing", "(4) Mov: CKMR+mixing")
)
readr::write_csv(Design, "tables/Design_10.02.2026_seasonal_mixing.csv")
