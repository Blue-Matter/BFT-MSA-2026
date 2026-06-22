
library(multiSA)

#### Calculate SD from a jittered run

fits <- readRDS(file = file.path("jitter", "jitter_reference1_06.03.2026.rds"))
i <- 5

fit <- fits[[i]]
dat <- get_MSAdata(fit)

parameters <- fit@obj$env$parList(fit@opt$par) # Estimated parameters in jittered run, no need to do optimization
random <- dat@Misc$random
map <- dat@Misc$map

fit_new <- fit_MSA(
  dat,
  parameters = parameters,
  random = random,
  map = map,
  run_model = FALSE,
  report = TRUE,
  do_sd = TRUE
)

# The starting parameters that generated this fit in the jitter run
fit_new@obj$par <- fit@obj$par
fit_new@opt <- fit@opt
saveRDS(fit_new, file = file.path("model_output", "newfit_reference1_06.03.2026.rds"))
