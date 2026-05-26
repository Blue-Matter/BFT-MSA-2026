
library(multiSA)

fit <- readRDS(file.path("model_output", "fit_reference_05.22.2026.rds"))
fitW <- readRDS(file.path("model_output", "fit_Wprior_05.22.2026.rds"))

# Single profile of WBFT R0
fit@report$R0_s
fitW@report$R0_s

tictoc::tic()
p <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = c(100, 150, 300, 500, 700, 1000, 1200, 1500, 2000, 3000),
  cores = 4
)
tictoc::toc()
saveRDS(p, file = file.path("profile", "WBFT_R0_05.22.2026.rds"))


#png("fit/figures/h_profile.png", height = 4, width = 5, res = 400, units = "in")
#par(mar = c(5, 4, 1, 1))
#plot(p, nlevels = 10, xlab = "EBFT steepness", ylab = "WBFT steepness", main = "Change in objective function")
#dev.off()
