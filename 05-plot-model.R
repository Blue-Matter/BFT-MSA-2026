
# plot
library(multiSA)
fit <- readRDS("model_output/fit_04.30.2026.rds")

dat <- get_MSAdata(fit)

dir_save <- file.path("figures", "fit")

# Aggregate fit to all indices
index_all <- lapply(1:dat@Dsurvey@ni, function(i) plot_index(fit, i = i, zoom = TRUE, figure = FALSE)) %>%
  bind_rows() %>%
  mutate(name = factor(name, dat@Dlabel@index))

g <- ggplot(index_all, aes(year, obs)) +
  geom_point(size = 0.75) +
  geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr)) +
  geom_line(aes(y = pred), colour = "red") +
  expand_limits(y = 0) +
  facet_wrap(vars(name), ncol = 3, scales = "free_y") +
  labs(x = "Year", y = "Index")
ggsave(file.path(dir_save, "index_fit.png"), g, width = 6, height = 8)

g <- ggplot(index_all, aes(year, log(obs/pred))) +
  geom_point(size = 0.75) +
  geom_line(linewidth = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(vars(name), ncol = 3, scales = "free_y") +
  labs(x = "Year", y = "log residual")
ggsave(file.path(dir_save, "index_resid.png"), g, width = 6, height = 8)

# Stock of origin fit
aa <- c("0-4", "5-8", "9+")
soo <- lapply(1:3, function(a) {
  plot_SC(fit, r = 2, a = a, prop = TRUE, figure = FALSE)
}) %>%
  bind_rows() %>%
  mutate(Age = .env$aa[.data$aa])

g <- soo %>%
  filter(stock == "WBFT") %>%
  filter(!is.na(pred), pred > 0) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(year, obs)) +
  geom_line(aes(y = pred), colour = 'red') +
  geom_line() +
  geom_point() +
  facet_grid(vars(region), vars(Age)) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_cartesian(xlim = c(1970, 2025)) +
  labs(x = "Year", y = "Proportion WBFT")
ggsave(file.path(dir_save, "SOO_fit.png"), g, height = 3, width = 6)

# Tags
png(file.path(dir_save, "tag_EBFT_fit.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 1)
dev.off()

png(file.path(dir_save, "tag_WBFT_fit.png"), height = 5, width = 6, res = 400, units = "in")
plot_tagmov(fit, s = 2)
dev.off()

# Mean length
mlen <- lapply(1:dat@Dfishery@nf, function(f) {
  lapply(1:dat@Dmodel@nr, function(r) {
    x <- plot_CAL(fit, f = f, r = r, do_mean = TRUE, figure = FALSE)
    if (is.null(x)) x <- data.frame()
    return(x)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

g <- mlen %>%
  filter(!is.na(pred)) %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet)) %>%
  mutate(region = factor(region, dat@Dlabel@region)) %>%
  arrange(year) %>%
  ggplot(aes(year, pred, fill = region, colour = region)) +
  geom_point(size = 1, shape = 21, colour = 'grey60', alpha = 0.75, aes(y = obs)) +
  geom_line() +
  facet_wrap(vars(fleet), ncol = 3) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "Mean length", region = NULL, colour = NULL)
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
  facet_wrap(vars(fleet), ncol = 3) +
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
  mutate(Fleet = factor(Fleet, dat@Dlabel@fleet))

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
png(file.path(dir_save, "recruitment.png"), height = 3, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
plot_R(fit, s = 1)
plot_R(fit, s = 2)
dev.off()

png(file.path(dir_save, "recdev.png"), height = 3, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
plot_Rdev(fit, s = 1)
title("EBFT")
plot_Rdev(fit, s = 2)
title("WBFT")
dev.off()

png(file.path(dir_save, "SRR.png"), height = 3, width = 8, units = "in", res = 400)
par(mfrow = c(1, 2))
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
plot_S(fit, by = "stock", facet_free = FALSE)
dev.off()

png(file.path(dir_save, "SSB_stock_compare.png"), height = 4, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "region", facet_free = FALSE)
dev.off()

png(file.path(dir_save, "SSB_stock_independent.png"), height = 4, width = 8, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot_S(fit, by = "region", facet_free = TRUE)
dev.off()

# Calculate regional exploitation rate
u <- local({
  CB_ymr <- apply(fit@report$CB_ymfrs, c(1, 2, 4), sum)
  B_ymr <- apply(fit@report$B_ymrs, 1:3, sum)
  Year = multiSA:::make_yearseason(dat@Dlabel@year, 4)
  U_yrs <- multiSA:::collapse_yearseason(CB_ymr/B_ymr) %>%
    structure(dimnames = list(Year = Year, Region = dat@Dlabel@region))
  reshape2::melt(U_yrs, value.name = "Ex")
}) %>%
  mutate(Season = 4 * (Year - floor(Year)) + 1)

#g <- ggplot(u, aes(Year, Ex)) +
#  facet_wrap(vars(Region)) +
#  geom_line(linewidth = 0.2) +
#  geom_point(alpha = 0.5, size = 0.75, aes(colour = factor(Season))) +
#  labs(y = "Seasonal Catch/Biomass", colour = "Season")
#ggsave(file.path(dir_save, "regional_exploitation.png"), g, height = 4, width = 5)

g <- u %>%
  mutate(Season = paste("Season", Season)) %>%
  ggplot(aes(floor(Year), Ex)) +
  facet_grid(vars(Season), vars(Region)) +
  geom_line() +
  #geom_point(alpha = 0.5, size = 0.75, aes(colour = factor(Season))) +
  labs(x = "Year", y = "Catch/Biomass")
ggsave(file.path(dir_save, "regional_exploitation.png"), g, height = 4, width = 6)

# Depletion
S_S0 <- local({
  S_ys <- apply(fit@report$S_yrs, c(1, 3), sum)
  S0 <- fit@report$SB0_s
  dep <- structure(t(S_ys)/S0, dimnames = list(Stock = dat@Dlabel@stock, Year = dat@Dlabel@year))

  reshape2::melt(dep, value.name = "dep")
})

g <- ggplot(S_S0, aes(Year, dep, colour = Stock)) +
  geom_line() +
  labs(y = expression(S/S[0]))
ggsave(file.path(dir_save, "depletion.png"), g, height = 3, width = 5)
