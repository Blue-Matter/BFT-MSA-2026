

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
# SSB_prior - logical, whether to use CKMR estimate of Western SSB in 2018
# SSB_sd - Numeric, decrease the SE of SSB prior
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
# (2a) With CKMR estimate, SD = 0.02
# (3) With CKMR + stock mixing in West Atlantic
Design <- data.frame(
  input_dir = "model_input/06.30.2026",
  annual = FALSE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(0, 0, 0, 1),
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = c(FALSE, TRUE, TRUE, TRUE),
  SSB_sd = c(0.18, 0.18, 0.02, 0.02),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = c(FALSE, FALSE, FALSE, TRUE),
  minI_SD = NA,
  Wareas = 2,
  Eareas = c(2, 2, 2, 3),
  output_name = paste0("seasonal_selprior", 1:4, "_08.19"),
  model_name = c("(1) NM: no CKMR", "(2) NM: CKMR SD = 0.18", "(2a) NM: CKMR SD = 0.02",
                 "(3) NM: CKMR+mixing")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal.csv")

# Same as Line 32 but with VAST indices
Design <- data.frame(
  input_dir = "model_input/06.30.2026",
  annual = FALSE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(0, 0, 1),
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = c(FALSE, TRUE, TRUE),
  SSB_sd = c(0.18, 0.02, 0.02),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = c(FALSE, FALSE, TRUE),
  minI_SD = NA,
  Wareas = 2,
  Eareas = c(2, 2, 3),
  output_name = paste0("seasonal_selprior", 1:3, "_08.19a"),
  model_name = c("(1) NM: no CKMR", "(2a) NM: CKMR SD = 0.02",
                 "(3) NM: CKMR+mixing")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_VAST.csv")

# Models with movement, first line with regular indices, second line with VAST indices
Design <- data.frame(
  input_dir = c("model_input/06.30.2026", "model_input/06.30.2026_VAST"),
  annual = FALSE,
  movement = TRUE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = 1,
  lambda_tag = 1,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = TRUE,
  SSB_sd = 0.02,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = TRUE,
  minI_SD = NA,
  Wareas = 2,
  Eareas = 3,
  output_name = paste0("seasonal_selprior", 1:2, "_08.19b"),
  model_name = "(4) Mov: CKMR+mixing"
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_movement.csv")


















Design <- data.frame(
  input_dir = "model_input/06.30.2026",
  annual = FALSE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(0, 0, 0, 1),
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = c(FALSE, TRUE, TRUE, TRUE),
  SSB_sd = c(0.18, 0.18, 0.02, 0.02),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = c(FALSE, FALSE, FALSE, TRUE),
  minI_SD = 0.1,
  Wareas = 2,
  Eareas = c(2, 2, 2, 3),
  output_name = paste0("seasonal_selprior", 1:4, "_09.16"),
  model_name = c("(1) NM: no CKMR", "(2) NM: CKMR SD = 0.18", "(2a) NM: CKMR SD = 0.02",
                 "(3) NM: CKMR+mixing")
)
readr::write_csv(Design, "tables/Design_09.16.2026_seasonal.csv")


# Models with movement, first line with regular indices, second line with VAST indices
Design <- data.frame(
  input_dir = c("model_input/06.30.2026", "model_input/06.30.2026_VAST"),
  annual = FALSE,
  movement = TRUE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = 1,
  lambda_tag = 1,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = TRUE,
  SSB_sd = 0.02,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = TRUE,
  minI_SD = c(0.1, NA),
  Wareas = 2,
  Eareas = 3,
  output_name = paste0("seasonal_selprior", 1:2, "_09.16"),
  model_name = "(4) Mov: CKMR+mixing"
)
readr::write_csv(Design, "tables/Design_09.16.2026_seasonal_movement.csv")


