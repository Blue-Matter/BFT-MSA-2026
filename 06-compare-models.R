
library(multiSA)
library(tidyverse)

source("99-functions-compare.R")

# Load design data frame of model fits
Design <- readr::read_csv("tables/Design_08.19.2026_seasonal.csv")[c(1, 2, 4), ] %>%
  mutate(model_name = c("(1) NM: no CKMR", "(2) NM: CKMR", "(3) NM: CKMR+SOO"))
dir_save <- "figures/fit/compare_08.19_seasonal"
table_suffix <- "08.19_seasonal"

Design <- rbind(
  readr::read_csv("tables/Design_08.19.2026_seasonal.csv")[1, ],
  readr::read_csv("tables/Design_08.19.2026_seasonal_CKMR02.csv")[c(1, 3), ]
) %>%
  mutate(model_name = c("(1) NM: no CKMR", "(2) NM: CKMR", "(3) NM: CKMR+SOO"))
dir_save <- "figures/fit/compare_08.19_seasonal_CKMR02"
table_suffix <- "08.19_seasonal_CKMR02"

Design <- readr::read_csv("tables/Design_08.19.2026_seasonal_VAST_CKMR02.csv")[c(1, 2, 4), ] %>%
  mutate(model_name = c("(1) NM: no CKMR", "(2) NM: CKMR", "(3) NM: CKMR+SOO"))
dir_save <- "figures/fit/compare_08.19_seasonal_VAST_CKMR02"
table_suffix <- "08.19_seasonal_VAST_CKMR02"

Design <- rbind(
  readr::read_csv("tables/Design_08.19.2026_seasonal_VAST_CKMR02.csv")[c(1, 2, 5), ],
  readr::read_csv("tables/Design_08.19.2026_seasonal_VAST_CKMR02_movement.csv")[1, ]
) %>%
  mutate(model_name = c("(1) NM: no CKMR", "(2) NM: CKMR", "(3) NM: CKMR+SOO", "(4) Mov: CKMR+SOO"))
dir_save <- "figures/fit/compare_08.19_seasonal_movement"
table_suffix <- "08.19_seasonal_movement"

Design <- readr::read_csv("tables/Design_08.19.2026_annual.csv")
dir_save <- "figures/fit/compare_08.19_annual"
table_suffix <- "08.19_annual"

if (!dir.exists(dir_save)) dir.create(dir_save)

# Load models
fits <- lapply(1:nrow(Design), function(i) {
  readRDS(file.path("model_output", paste0("fit_", Design$output_name[i], ".rds")))
})
dat <- get_MSAdata(fits[[1]])



# Diagnostic stuff
if (FALSE) {
  # Adjacency matrix

  library(igraph)

  H <- fits[[1]]@SD$env$hessian
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
}



# Compare SSB in spawning season ----
# Spawning biomass - common y-axis between stocks
g <- plot_SSB(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_SSB.png"), g, height = 5, width = 6)

# Spawning biomass - separate axes by stock
g <- plot_SSB(fits, Design$model_name, scales = "free_y")
ggsave(file.path(dir_save, "compare_SSB2.png"), g, height = 5, width = 6)

# Recruitment ----
g <- plot_rec(fits, Design$model_name)
ggsave(file.path(dir_save, "compare_rec.png"), g, height = 4, width = 6)

g <- plot_recdev(fits, Design$model_name) +
  geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr))
#g <- plot_recdev(fits, Design$model_name)
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
g <- plot_fit_soo3(fits, Design$model_name) +
  guides(colour = guide_legend(ncol = 2))
ggsave(file.path(dir_save, "compare_SOO3_fit.png"), g, height = 5, width = 8)

# CAL ----
g <- plot_fit_CAL_agg(fits, Design$model_name) +
  guides(colour = guide_legend(ncol = 2))
ggsave(file.path(dir_save, "compare_CAL_agg_fit.png"), g, height = 8, width = 6)

#summarise(g@data, CAL_mode = lmid[which.max(obs)], .by = fleet) %>%
#  mutate(LFS = fits[[2]]@report$selconv_pf[1, 1:18],
#         init = fits[[2]]@obj$report(fits[[2]]@obj$par)$selconv_pf[1, 1:18])

# Mean length ----
g <- plot_mlen(fits, Design$model_name) +
  guides(colour = guide_legend(ncol = 2))
ggsave(file.path(dir_save, "compare_mlen_fit.png"), g, height = 8, width = 6)


# Only for models with stock-specific selectivity
if (any(Design$est_stocksel)) {
  j <- which(Design$est_stocksel)
  j <- 4

  # Plot F by stock ----
  g <- fits[[j]]@report$F_ymafrs %>%
    apply(c(1, 2, 4, 5, 6), max) %>%
    structure(
      dimnames = list(year = dat@Dlabel@year,
                      season = seq(1, dat@Dmodel@nm),
                      fleet = dat@Dlabel@fleet,
                      region = dat@Dlabel@region,
                      stock = dat@Dlabel@stock)
    ) %>%
    reshape2::melt(id.vars = "F") %>%
    filter(region == "WATL") %>%
    filter(sum(value) > 0, .by = fleet) %>%
    summarise(value = sum(value), .by = c(year, fleet, stock, region)) %>%
    #filter(diff(value) != 0, .by = c(year, fleet, region)) %>%
    ggplot(aes(year, value, linetype = stock)) +
    geom_line() +
    facet_wrap(vars(fleet), scales = "free_y") +
    labs(x = "Year", y = "Fishing mortality (per year)", linetype = NULL, title = "WATL") +
    theme(legend.position = "bottom")
  ggsave(file.path(dir_save, "F_stock_WATL.png"), g, height = 5, width = 6)

}



# Only for models with movement
if (any(Design$movement)) {
  j <- which(Design$movement)
  j <- 4

  # Plot tags ----
  g <- plot_tags(fits[j], Design$model_name[j], stock = "EBFT") +
    labs(title = "EBFT, Age 0-4")
  ggsave(file.path(dir_save, "compare_tag_EBFT_a1.png"), g, height = 6, width = 6)

  g <- plot_tags(fits[j], Design$model_name[j], stock = "EBFT", ac = 2) +
    labs(title = "EBFT, Age 5-8")
  ggsave(file.path(dir_save, "compare_tag_EBFT_a2.png"), g, height = 6, width = 6)

  g <- plot_tags(fits[j], Design$model_name[j], stock = "EBFT", ac = 3) +
    labs(title = "EBFT, Age 9+")
  ggsave(file.path(dir_save, "compare_tag_EBFT_a3.png"), g, height = 6, width = 6)

  g <- plot_tags(fits[j], Design$model_name[j], stock = "WBFT", ac = 3) +
    labs(title = "WBFT, Age 9+")
  ggsave(file.path(dir_save, "compare_tag_WBFT_a3.png"), g, height = 6, width = 6)

  ## Plot movement ----
  g <- plot_mov_gg(fits[j], Design$model_name[j])
  g <- plot_mov_gg(fits[j], Design$model_name[j], stock = "WBFT")

  png(file.path(dir_save, "mov_EBFT.png"), res = 400, width = 8, height = 6, units = "in")
  par(mar = c(5, 4, 1, 1))
  plot_mov(fits[[j]], s = 1)
  dev.off()

  png(file.path(dir_save, "mov_WBFT.png"), res = 400, width = 8, height = 6, units = "in")
  par(mar = c(5, 4, 1, 1))
  plot_mov(fits[[j]], s = 2)
  dev.off()

  # Plot mature biomass by area & season ----
  g <- plot_SB_season(fits[[j]]) +
    labs(title = Design$model_name[j])
  ggsave(file.path(dir_save, "SB_area_season.png"), g, height = 5, width = 6)

  # Plot mature biomass by area & season ----
  g@facet$params$free$y <- TRUE
  ggsave(file.path(dir_save, "SB_area_season2.png"), g, height = 5, width = 6)

}

# Compare SSB with M3 and SS3 ----
.g <- plot_SSB(fits, Design$model_name)
M3 <- readr::read_csv("tables/M3_SSB.csv") %>%
  mutate(S = 1e-3 * S)
SS3 <- readr::read_csv("tables/SS3_SSB.csv") %>%
  mutate(stock = ifelse(area == "WATL", "WBFT", "EBFT"))

SSB <- .g@data %>%
  mutate(model = paste("MSA:", model)) %>%
  summarise(S = sum(S), .by = c(year, stock, model)) %>%
  rbind(dplyr::select(M3, year, stock, S, model)) %>%
  rbind(dplyr::select(SS3, year, stock, S, model)) %>%
  mutate(software = ifelse(grepl("M3", model), "M3",
                           ifelse(grepl("SS3", model), "SS3", "MSA"))) %>%
  arrange(year)

g <- ggplot(SSB, aes(year, S, colour = model, group = model)) +
  geom_line() +
  facet_grid(vars(stock), vars(software), scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Spawning biomass", colour = NULL) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 4))
ggsave(file.path(dir_save, "compare_SSB_M3_SS3.png"), g, height = 5, width = 6)

