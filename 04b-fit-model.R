

library(multiSA)
library(tidyverse)
library(tictoc)
library(parallel)

source("99-functions-fit.R")

# Load design data frame
#Design <- readr::read_csv("tables/Design_08.19.2026_seasonal.csv")
#Design <- readr::read_csv("tables/Design_08.19.2026_seasonal_CKMR02.csv")
#Design <- rbind(
#  readr::read_csv("tables/Design_08.19.2026_seasonal_VAST_CKMR02.csv"),
#  readr::read_csv("tables/Design_08.19.2026_seasonal_VAST_CKMR02_movement.csv")
#)
Design <- readr::read_csv("tables/Design_08.19.2026_annual.csv")

# Fit all models in parallel or in a loop ----
do_parallel <- TRUE

if (do_parallel) {

  cl <- parallel::makeCluster(min(nrow(Design), 4)) # Need to limit the number of cores due to memory constraints
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

