
library(multiSA)
library(tidyverse)
library(randtests)

#fit <- readRDS("model_output/fit_reference_05.22.2026.rds")
#dir_save <- file.path("figures", "fit", "05.22")

fit <- readRDS("model_output/fit_Wprior_05.22.2026.rds")
dir_save <- file.path("figures", "fit", "05.22_Wprior")
if (!dir.exists(dir_save)) dir.create(dir_save)

dat <- get_MSAdata(fit)

# Aggregate fit to all indices
index_all <- lapply(1:dat@Dsurvey@ni, function(i) plot_index(fit, i = i, zoom = TRUE, figure = FALSE)) %>%
  bind_rows() %>%
  mutate(resid = log(obs/pred)) %>%
  mutate(name = factor(name, dat@Dlabel@index))

run_test <- data.frame(
  name = dat@Dlabel@index,
  pass = sapply(dat@Dlabel@index, function(i) {
    p.value <- filter(index_all, name == i) %>%
      pull(resid) %>%
      randtests::runs.test(alternative = "left.sided", threshold = 0, plot = FALSE) %>%
      getElement("p.value")
    ifelse(p.value > 0.05, "Pass", "Fail")
  })
) %>%
  mutate(col = ifelse(pass == "Pass", "black", "red")) %>%
  mutate(name = factor(name, dat@Dlabel@index))

g <- ggplot(index_all, aes(year, obs)) +
  geom_point(size = 0.75) +
  geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr)) +
  geom_line(aes(y = pred), colour = "red") +
  expand_limits(y = 0) +
  facet_wrap(vars(name), ncol = 3, scales = "free_y") +
  labs(x = "Year", y = "Index")
ggsave(file.path(dir_save, "index_fit.png"), g, width = 6, height = 8)

g <- ggplot(index_all, aes(year, resid)) +
  geom_label(data = run_test, x = -Inf, y = Inf, hjust = "inward", vjust = "inward", aes(colour = col, label = pass)) +
  geom_point(size = 0.75) +
  geom_line(linewidth = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(vars(name), ncol = 3, scales = "free_y") +
  scale_colour_identity() +
  labs(x = "Year", y = "log residual")
ggsave(file.path(dir_save, "index_resid.png"), g, width = 6, height = 8)

# Stock of origin fit
aa <- c("0-4", "5-8", "9+")
soo <- lapply(1:3, function(a) {
  lapply(2:3, function(r) {
    lapply(1:2, function(ff) {
      plot_SC(fit, r = r, ff = ff, a = a, prop = TRUE, figure = FALSE)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(Age = .env$aa[.data$aa],
         region = factor(region, c("WATL", "EATL")),
         Fleet = ifelse(ff == 1, "Otolith", "Genetic"))

g <- soo %>%
  filter(stock == "WBFT") %>%
  filter(!is.na(pred), pred > 0) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(year, obs, colour = Fleet, fill = Fleet, shape = Fleet)) +
  geom_line(linewidth = 0.1) +
  geom_line(aes(y = pred), linewidth = 0.25, colour = "black") +
  #geom_line(aes(y = pred), colour = 'red') +
  geom_point() +
  facet_grid(vars(region), vars(Age)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_shape_manual(values = c(1, 8)) +
  #coord_cartesian(xlim = c(1970, 2025)) +
  labs(x = "Year", y = "Proportion WBFT", fill = NULL, shape = NULL, colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(dir_save, "SOO_fit.png"), g, height = 5, width = 8)


SOO_resid <- structure(
  residuals(fit, vars = "SC_ymafrs")$SC_ymafrs,
  dimnames = list(
    Year = dat@Dlabel@year,
    Season = 1:dat@Dmodel@nm,
    Age = aa,
    Fleet = c("Otolith", "Genetic"),
    Region = dat@Dlabel@region,
    Stock = dat@Dlabel@stock
  )
) %>%
  reshape2::melt()
g <- SOO_resid %>%
  filter(Region %in% c("WATL", "EATL")) %>%
  #filter(!is.na(value)) %>%
  filter(Stock == "WBFT") %>%
  mutate(Age = paste("Age", Age),
         Fleet = factor(Fleet, c("Genetic", "Otolith")),
         year = Year + 0.25 * (Season - 1)) %>%
  ggplot(aes(year, value, shape = Fleet, colour = Fleet)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(1, 8)) +
  facet_grid(vars(Region), vars(Age)) +
  labs(x = "Year", y = "SOO residual") +
  theme(legend.position = "bottom")
ggsave(file.path(dir_save, "SOO_residual.png"), g, height = 5, width = 8)


# Tags
png(file.path(dir_save, "tag_EBFT_fit_a1.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 1, aa = 1)
dev.off()

png(file.path(dir_save, "tag_EBFT_fit_a2.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 1, aa = 2)
dev.off()

png(file.path(dir_save, "tag_EBFT_fit_a3.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 1, aa = 3)
dev.off()

png(file.path(dir_save, "tag_WBFT_fit.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 2, aa = 3)
dev.off()

# Mean length
panel_factor <- outer(dat@Dlabel@fleet, dat@Dlabel@region, function(i, j) paste(i, j, sep = ": ")) %>%
  t() %>%
  as.character()

mlen <- lapply(1:dat@Dfishery@nf, function(f) {
  lapply(1:dat@Dmodel@nr, function(r) {
    x <- plot_CAL(fit, f = f, r = r, do_mean = TRUE, figure = FALSE)
    if (is.null(x)) x <- data.frame()
    return(x)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet)) %>%
  mutate(region = factor(region, dat@Dlabel@region)) %>%
  arrange(year) %>%
  mutate(panel = paste(fleet, region, sep = ": ") %>% factor(panel_factor))

g <- mlen %>%
  filter(!is.na(pred)) %>%
  ggplot(aes(year, pred)) +
  geom_point(alpha = 0.5, shape = 1, aes(y = obs)) +
  geom_line(data = filter(mlen, !is.na(obs)), linewidth = 0.1, aes(y = obs)) +
  geom_line(colour = "red") +
  facet_wrap(vars(panel), ncol = 4) +
  #expand_limits(y = 0) +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "Mean length", fill = NULL, colour = NULL)
ggsave(file.path(dir_save, "mlen.png"), g, height = 8, width = 6)


# Aggregate length comp
CAL <- lapply(1:dat@Dfishery@nf, function(f) {
  lapply(1:dat@Dmodel@nr, function(r) {
    x <- plot_CAL(fit, f = f, r = r, do_mean = FALSE, figure = FALSE)
    if (is.null(x)) x <- data.frame()
    return(x)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

CAL_agg <- CAL %>%
  mutate(  # for each time step, fleet, and region, re-do N
    obs2 = N * obs/sum(obs, na.rm = TRUE),
    pred2 = N * pred/sum(pred, na.rm = TRUE),
    .by = c(year, fleet, region)
  ) %>%
  summarise(
    obs = sum(obs2, na.rm = TRUE),
    pred = sum(pred2, na.rm = TRUE),
    .by = c(fleet, lmid)
  ) %>%
  mutate(
    obs = obs/sum(obs),
    pred = pred/sum(pred),
    .by = fleet
  )

g <- CAL_agg %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet)) %>%
  ggplot(aes(lmid, obs)) +
  geom_area(fill = "grey80", linewidth = 0.1, colour = "black") +
  geom_point() +
  geom_line(aes(y = pred), linewidth = 1.25, colour = "mediumseagreen") +
  facet_wrap(vars(fleet), ncol = 3, scales = "free_y") +
  labs(x = "Length Bin", y = "Proportion")
ggsave(file.path(dir_save, "CAL_agg_fit.png"), g, height = 8, width = 6)

# Residuals CAL
CAL_resid <- structure(residuals(fit, vars = "CALobs_ymlfr")$CALobs_ymlfr,
  dimnames = list(
    Year = dat@Dlabel@year,
    Season = dat@Dlabel@season,
    Length = dat@Dmodel@lmid,
    Fleet = dat@Dlabel@fleet,
    Region = dat@Dlabel@region
  )) %>%
  reshape2::melt() %>%
  mutate(Fleet = factor(Fleet, dat@Dlabel@fleet)) %>%
  filter(!is.na(value)) %>%
  mutate(value = pmin(value, 5) %>% pmax(-5))

g <- ggplot(CAL_resid) +
  geom_histogram(fill = "grey60", colour = "black", linewidth = 0.1,
                 aes(value, after_stat(density))) +
  facet_wrap(vars(Fleet), ncol = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(x = "CAL residuals")
ggsave(file.path(dir_save, "CAL_resid_hist.png"), g, height = 8, width = 6)

# Plot selectivity
self <- lapply(1:dat@Dfishery@nf, function(f) {
  plot_self(fit, f = f, figure = FALSE)
}) %>%
  bind_rows() %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet))

g <- ggplot(self, aes(length, sel)) +
  geom_line() +
  facet_wrap(vars(fleet), ncol = 3) +
  labs(x = "Length", y = "Selectivity")
ggsave(file.path(dir_save, "sel_len.png"), g, height = 8, width = 6)

self_age <- lapply(1:dat@Dfishery@nf, function(f) {
  plot_self(fit, f = f, type = "age", figure = FALSE)
}) %>%
  bind_rows() %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet))

g <- ggplot(self_age, aes(age, sel)) +
  geom_line() +
  facet_wrap(vars(fleet), ncol = 3) +
  labs(x = "Age", y = "Selectivity")
ggsave(file.path(dir_save, "sel_age.png"), g, height = 8, width = 6)

seli <- lapply(1:dat@Dsurvey@ni, function(i) {
  plot_seli(fit, i = i, figure = FALSE)
}) %>%
  bind_rows() %>%
  mutate(name = factor(name, dat@Dlabel@index))

g <- seli %>%
  filter(!is.na(length)) %>%
  ggplot(aes(length, sel)) +
  geom_line() +
  facet_wrap(vars(name), ncol = 4) +
  labs(x = "Length", y = "Selectivity")
ggsave(file.path(dir_save, "selindex_length.png"), g, height = 5, width = 6)

g <- seli %>%
  filter(is.na(length)) %>%
  ggplot(aes(age, sel)) +
  geom_line() +
  facet_wrap(vars(name)) +
  labs(x = "Age", y = "Selectivity")


# Recruitment
png(file.path(dir_save, "recruitment.png"), height = 4, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
plot_R(fit, s = 1)
plot_R(fit, s = 2)
dev.off()

png(file.path(dir_save, "recdev.png"), height = 4, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
plot_Rdev(fit, s = 1)
title("EBFT")
plot_Rdev(fit, s = 2)
title("WBFT")
dev.off()

png(file.path(dir_save, "SRR.png"), height = 4, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
plot_SRR(fit, s = 1)
title("EBFT")
plot_SRR(fit, s = 2)
title("WBFT")
dev.off()

# Recdist
png(file.path(dir_save, "recdist.png"), height = 4, width = 4, units = "in", res = 400)
par(mfrow = c(1, 1), mar = c(5, 4, 1, 1))
plot_recdist(fit)
dev.off()

# Movement
png(file.path(dir_save, "mov_EBFT.png"), height = 6, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_mov(fit, s = 1)
dev.off()

png(file.path(dir_save, "mov_WBFT.png"), height = 6, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_mov(fit, s = 2)
dev.off()

# SSB
png(file.path(dir_save, "SSB_area.png"), height = 6, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "stock", facet_free = FALSE, ylab = "Spawning stock biomass")
dev.off()

png(file.path(dir_save, "SSB_area_prop.png"), height = 6, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "stock", facet_free = FALSE, prop = TRUE, ylab = "Spawning stock biomass")
dev.off()

png(file.path(dir_save, "SSB_stock_compare.png"), height = 4, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "region", facet_free = FALSE, ylab = "Spawning stock biomass")
dev.off()

png(file.path(dir_save, "SSB_stock_independent.png"), height = 4, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "region", facet_free = TRUE, ylab = "Spawning stock biomass")
dev.off()

# Spawning biomass by season
SB_season <- local({
  N_ymars <- fit@report$N_ymars[1:dat@Dmodel@ny, , , , ]
  mat_ymars <- array(
    dat@Dstock@mat_yas,
    c(dat@Dmodel@ny, dat@Dmodel@na, dat@Dmodel@ns, dat@Dmodel@nm, dat@Dmodel@nr)
  ) %>%
    aperm(c(1, 4, 2, 5, 3))
  fec_ymars <- array(
    dat@Dstock@swt_ymas,
    c(dat@Dmodel@ny, dat@Dmodel@nm, dat@Dmodel@na, dat@Dmodel@ns, dat@Dmodel@nr)
  ) %>%
    aperm(c(1, 2, 3, 5, 4))

  S_ymrs <- apply(N_ymars * mat_ymars * fec_ymars, c(1, 2, 4, 5), sum) %>%
    structure(dimnames = list(
      year = dat@Dlabel@year,
      season = dat@Dlabel@season,
      region = dat@Dlabel@region,
      stock = dat@Dlabel@stock
    ))
  reshape2::melt(S_ymrs, value.name = "S")
})
g <- SB_season %>%
  ggplot(aes(year, S, fill = stock)) +
  geom_col(width = 1) +
  facet_grid(vars(region), vars(season)) +
  labs(x = "Year", y = "Mature biomass", fill = NULL) +
  scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dir_save, "SSB_area_season.png"), g, height = 5, width = 6)



# Biomass by season
png(file.path(dir_save, "B_area.png"), height = 6, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_B(fit, by = "stock", facet_free = FALSE)
dev.off()

g <- plot_B(fit, by = "stock", facet_free = FALSE, figure = FALSE) %>%
  mutate(season = 4 * (year - floor(year)) + 1) %>%
  mutate(season = paste("Season", season)) %>%
  ggplot(aes(floor(year), B, fill = stock)) +
  geom_col(width = 1) +
  facet_grid(vars(region), vars(season)) +
  labs(x = "Year", y = "Total biomass", fill = NULL) +
  scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dir_save, "B_area_season.png"), g, height = 5, width = 6)

# Calculate regional exploitation rate
u <- local({
  CB_ymrs <- apply(fit@report$CB_ymfrs, c(1, 2, 4, 5), sum)
  B_ymrs <- fit@report$B_ymrs
  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
  U_ymrs <- CB_ymrs/B_ymrs
  U_ymrs[CB_ymrs < 1e-8] <- 0
  U_yrs <- multiSA:::collapse_yearseason(U_ymrs) %>%
    structure(dimnames = list(Year = Year, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock))
  reshape2::melt(U_yrs, value.name = "Ex")
}) %>%
  mutate(Season = 4 * (Year - floor(Year)) + 1)

u_eq <- local({
  CB_mrs <- apply(fit@report$initCB_mfrs, c(1, 3, 4), sum)
  N_mars <- sapply2(1:dat@Dmodel@ns, function(s) fit@report$initNPR_mars[, , , s] * fit@report$initReq_s[s])
  B_mrs <- sapply2(1:dat@Dmodel@nr, function(r) N_mars[, , r, ] * dat@Dstock@swt_ymas[1, , , ]) %>%
    apply(c(1, 4, 3), sum)

  U_mrs <- CB_mrs/B_mrs
  U_mrs[CB_mrs < 1e-8] <- 0

  structure(U_mrs, dimnames = list(Season = dat@Dlabel@season, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock)) %>%
    reshape2::melt(value.name = "Ex")
}) %>%
  mutate(Year = min(u$Year) - 1)

g <- u %>%
  mutate(Season = paste("Season", Season)) %>%
  ggplot(aes(floor(Year), Ex, colour = Stock)) +
  facet_grid(vars(Region), vars(Season), scales = "free") +
  geom_line() +
  geom_point(data = u_eq) +
  #geom_point(alpha = 0.5, size = 0.75, aes(colour = factor(Season))) +
  labs(x = "Year", y = "Seasonal Catch/Biomass") +
  theme(legend.position = "bottom")
ggsave(file.path(dir_save, "regional_exploitation.png"), g, height = 5, width = 6)

# Apical fishing mortality by stock
plot_Fstock(fit, s = 1:2, 'season')

# Depletion
S_S0 <- local({
  S_ys <- apply(fit@report$S_yrs, c(1, 3), sum)
  S0 <- fit@report$SB0_s
  dep <- structure(t(S_ys)/S0, dimnames = list(Stock = dat@Dlabel@stock, Year = dat@Dlabel@year))

  reshape2::melt(dep, value.name = "dep")
})

g <- ggplot(S_S0, aes(Year, dep, linetype = Stock)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1.5)) +
  expand_limits(y = 0) +
  labs(y = expression(S/S[0]))
ggsave(file.path(dir_save, "depletion.png"), g, height = 3, width = 5)

