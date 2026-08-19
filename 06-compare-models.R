
library(multiSA)
library(tidyverse)

source("99-functions-compare.R")

Design <- readr::read_csv("tables/Design_08.19.2026_annual.csv")
#Design <- readr::read_csv("tables/Design_08.04.2026.csv")[4:7, ]
Design <- readr::read_csv("tables/Design_08.12.2026.csv")[1:2, ]


fits <- lapply(1:nrow(Design), function(i) {
  readRDS(file.path("model_output", paste0("fit_", Design$output_name[i], ".rds")))
})
dat <- get_MSAdata(fits[[1]])

dir_save <- "figures/fit/compare_08.19_annual"
table_suffix <- "08.12_annual"

# Adjacency matrix
library(igraph)

H <- fits[[3]]@SD$env$hessian
A <- (abs(H) > 0)  # adjacency matrix

g <- graph_from_adjacency_matrix(A, mode = "undirected")
comp <- components(g)
plot(g, vertex.color = comp$membership)


# Convergence ----
lapply(fits, function(i) i@SD$pdHess)

# Likelihoods ----
like <- lapply(fits, multiSA:::get_likelihood_components) %>%
  bind_rows() %>%
  as.matrix() %>%
  `rownames<-`(Design$model_name); like
write.csv(as.data.frame(like), paste0("tables/like_", table_suffix, ".csv"))

like %>%
  apply(2, function(x) x - max(x)) %>%
  t() %>%
  signif(4)

# Likelihood of CAL by fleet
sapply(fits, function(i) {
  apply(i@report$loglike_CAL_ymfr, 3, sum)
})
apply(dat@Dfishery@CALN_ymfr, 3, sum) # Sample size

# Correlation ----
corr <- get_corr(fits[[1]])
condition_number(fits[[1]])

Francis_weights(fits[[1]]) %>%
  mutate(w = round(w, 2))

# Get values of R0
var_name <- "R0_s"
R0_df <- sapply(fits, function(i) i@report[[var_name]]) %>%
  structure(dimnames = list(Stock = dat@Dlabel@stock, Model = Design$model_name)) %>%
  t() %>%
  as.data.frame() %>%
  round() %>%
  mutate(ratio = round(EBFT/WBFT, 2))
write.csv(R0_df, file = paste0("tables/compare_R0_", table_suffix, ".csv"))

# Compare SSB in spawning season ----
# Spawning biomass - common y-axis between stocks
g <- plot_SSB(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_SSB.png"), g, height = 5, width = 7)

# Spawning biomass - separate axes by stock
g <- plot_SSB(fits, Design$model_name, scales = "free_y")
ggsave(file.path(dir_save, "compare_SSB2.png"), g, height = 5, width = 7)

# Recruitment ----
g <- plot_rec(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_rec.png"), g, height = 4, width = 6)

g <- plot_recdev(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_recdev.png"), g, height = 4, width = 6)

# Selectivity ----
g <- plot_sel(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_sel.png"), g, height = 8, width = 6)

g <- plot_sel_index(fits, Design$model_name, type = "length")
ggsave(file.path(dir_save, "compare_sel_cpue.png"), g, height = 8, width = 6)

g <- plot_sel_index(fits, Design$model_name, type = "age")
ggsave(file.path(dir_save, "compare_sel_index.png"), g, height = 5, width = 6)

# Fit to indices ----
# CPUE
g <- plot_fit_index(fits, Design$model_name, "cpue")
ggsave(file.path(dir_save, "compare_CPUE_fit.png"), g, width = 6, height = 8)

# Fishery-independent indices
g <- plot_fit_index(fits, Design$model_name, "fi")
ggsave(file.path(dir_save, "compare_index_fit.png"), g, width = 6, height = 6)

# SOO ----
g <- plot_fit_soo3(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_SOO3_fit.png"), g, height = 5, width = 8)

# CAL ----
g <- plot_fit_CAL_agg(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_CAL_agg_fit.png"), g, height = 8, width = 6)

summarise(g@data, CAL_mode = lmid[which.max(obs)], .by = fleet) %>%
  mutate(LFS = fits[[2]]@report$selconv_pf[1, 1:18],
         init = fits[[2]]@obj$report(fits[[2]]@obj$par)$selconv_pf[1, 1:18])

# Mean length ----
g <- plot_mlen(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_mlen_fit.png"), g, height = 8, width = 6)
