

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
  SSB_sd = 0.18,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  Wareas = 2,
  Eareas = c(2, 2, 3, 3),
  output_name = paste0("seasonal_selprior", 1:4, "_08.19"),
  model_name = c("1. NM, No CKMR", "2. NM, CKMR", "3. NM, CKMR + WATL mixing", "4. NM, (3) + SOO")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal.csv")

# Seasonal, no movement, tighter CKMR cv
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
  SSB_prior = TRUE,
  SSB_sd = 0.02,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  Wareas = 2,
  Eareas = c(2, 3, 3),
  output_name = paste0("seasonal_selprior", 2:4, "_08.19a"),
  model_name = c("2. NM, CKMR", "3. NM, CKMR + WATL mixing", "4. NM, (3) + SOO")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_CKMR02.csv")

Design <- data.frame(
  input_dir = "model_input/06.30.2026_VAST",
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
  SSB_sd = 0.18,
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  Wareas = 2,
  Eareas = c(2, 2, 3, 3),
  output_name = paste0("seasonal_selprior", 1:4, "_08.19c"),
  model_name = c("1. NM, No CKMR", "2. NM, CKMR", "3. NM, CKMR + WATL mixing", "4. NM, (3) + SOO")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_VAST.csv")

Design <- data.frame(
  input_dir = "model_input/06.30.2026_VAST",
  annual = FALSE,
  movement = FALSE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(0, 0, 0, 1, 1),
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = c(FALSE, TRUE, TRUE, TRUE, TRUE),
  SSB_sd = c(0.18, 0.02, 0.02, 0.02, 0.02),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  Wareas = 2,
  Eareas = c(2, 2, 3, 3, 3),
  output_name = paste0("seasonal_selprior", 1:5, "_08.19d"),
  model_name = c("1. NM, No CKMR", "2. NM, CKMR", "3. NM, CKMR + WATL mixing", "4. NM, (3) + SOO", "5. NM, (3) + SOO stocksel")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_VAST_CKMR02.csv")

Design <- data.frame(
  input_dir = "model_input/06.30.2026_VAST",
  annual = FALSE,
  movement = TRUE,
  rec_devvector = FALSE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(1, 1),
  lambda_tag = c(1, 1),
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = TRUE,
  SSB_sd = 0.02,
  spat_prior = c(FALSE, TRUE),
  sel_prior = TRUE,
  fix_sel = FALSE,
  est_stocksel = TRUE,
  Wareas = 2,
  Eareas = 3,
  output_name = paste0("seasonal_selprior", 6:7, "_08.19d"),
  model_name = c("6. Mov, SOO + Tag", "7. Mov, (6) + spat. prior")
)
readr::write_csv(Design, "tables/Design_08.19.2026_seasonal_VAST_CKMR02_movement.csv")


# Annual model
Design <- data.frame(
  input_dir = c("model_input/06.30.2026_annual_2area", "model_input/06.30.2026_annual_4area",
                "model_input/06.30.2026_annual_2area", "model_input/06.30.2026_annual_4area"),
  annual = TRUE,
  movement = FALSE,
  rec_devvector = TRUE,
  initC_scalar = 0.5,
  lambda_CAL = 1,
  lambda_SC = c(0, 0, 0, 0),
  lambda_tag = 0,
  SC_set = 3,
  SC_subset = "all",
  SSB_prior = c(FALSE, FALSE, TRUE, TRUE),
  SSB_sd = c(0.18, 0.18, 0.02, 0.02),
  spat_prior = FALSE,
  sel_prior = TRUE,
  fix_sel = FALSE,
  Wareas = c(1, 2, 1, 2),
  Eareas = c(1, 2, 1, 2),
  output_name = c("annual_selprior1_2area_08.19", "annual_selprior1_4area_08.19",
                  "annual_selprior2_2area_08.19", "annual_selprior2_4area_08.19"),
  model_name = c("Annual 2A", "Annual 4A", "Annual 2A+CKMR", "Annual, 4A+CKMR")
)
readr::write_csv(Design, "tables/Design_08.19.2026_annual.csv")

