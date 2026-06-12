
library(multiSA)
library(tidyverse)

Design <- readr::read_csv("tables/Design_06.03.2026.csv")[1:4, ]

#Design$output_name[c(2, 4)] <- paste0(Design$output_name[c(2, 4)], "_CV05")

fits <- lapply(1:nrow(Design), function(i) {
  if (i == 1) {
    readRDS(file.path("model_output", paste0("newfit_", Design$output_name[i], ".rds")))
  } else {
    readRDS(file.path("model_output", paste0("fit_", Design$output_name[i], ".rds")))
  }
})


# Convergence
lapply(fits, function(i) i@SD$pdHess)

# Likelihoods
like <- lapply(fits, multiSA:::get_likelihood_components) %>%
  bind_rows() %>%
  as.matrix() %>%
  `rownames<-`(Design$model_name)

like %>%
  apply(2, function(x) x - max(x)) %>%
  t() %>%
  signif(4)

# Get values of parameters
dat <- get_MSAdata(fits[[1]])

var_name <- "R0_s"
R0_df <- sapply(fits, function(i) i@report[[var_name]]) %>%
  structure(dimnames = list(Stock = dat@Dlabel@stock, Model = Design$model_name)) %>%
  t() %>%
  as.data.frame() %>%
  round() %>%
  mutate(ratio = round(EBFT/WBFT, 2))
write.csv(R0_df, file = "tables/compare_R0.csv")

# Compare SSB
SSB <- lapply(1:nrow(Design), function(i) {
  plot_S(fits[[i]], s = 1:2, figure = FALSE) %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, Design$model_name))

prior <- data.frame(
  model = Design$model_name[c(2, 4)],
  stock = "WBFT",
  year = 2018,
  S = 22000,
  lwr = exp(log(22000) - 1.96 * 0.18),
  upr = exp(log(22000) + 1.96 * 0.18)
) %>%
  mutate(model = factor(model, Design$model_name))

g <- SSB %>%
  ggplot(aes(year, S)) +
  facet_grid(vars(stock), vars(model)) +
  geom_col(width = 1, aes(fill = region)) +
  geom_pointrange(data = prior, size = 0.25, aes(ymin = lwr, ymax = upr)) +
  labs(x = "Year", y = "Spawning stock biomass (season 2)", fill = NULL) +
  scale_fill_manual(values = multiSA:::make_color(4, "region")) +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_SSB.png", g, height = 5, width = 7)

g <- SSB %>%
  ggplot(aes(year, S)) +
  facet_grid(vars(stock), vars(model), scales = "free_y") +
  geom_col(width = 1, aes(fill = region)) +
  geom_pointrange(data = prior, size = 0.25, aes(ymin = lwr, ymax = upr)) +
  labs(x = "Year", y = "Spawning stock biomass (season 2)", fill = NULL) +
  scale_fill_manual(values = multiSA:::make_color(4, "region")) +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_SSB2.png", g, height = 5, width = 7)





# Aggregate fit to all indices
index_all <- lapply(1:length(fits), function(i) {
  dat <- get_MSAdata(fits[[i]])

  index_all <- lapply(1:dat@Dsurvey@ni, function(ii) plot_index(fits[[i]], i = ii, zoom = TRUE, figure = FALSE)) %>%
    bind_rows() %>%
    mutate(resid = log(obs/pred)) %>%
    mutate(name = factor(name, dat@Dlabel@index)) %>%
    mutate(model = Design$model_name[i])

  return(index_all)
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, Design$model_name))

g <- ggplot(index_all, aes(year, obs)) +
  geom_point(size = 0.5) +
  geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr)) +
  geom_line(aes(y = pred, colour = model)) +
  expand_limits(y = 0) +
  facet_wrap(vars(name), ncol = 3, scales = "free_y") +
  labs(x = "Year", y = "Index", colour = "Model") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
ggsave("figures/fit/compare_index_fit.png", g, width = 6, height = 8)

# SOO
aa <- c("0-4", "5-8", "9+")
soo <- lapply(1:length(fits), function(i) {
  dat <- get_MSAdata(fits[[i]])

  lapply(1:3, function(a) {
    lapply(2:3, function(r) {
      lapply(1:2, function(ff) {
        plot_SC(fits[[i]], r = r, ff = ff, aa = a, prop = TRUE, figure = FALSE)
      }) %>%
        bind_rows()
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(Age = .env$aa[.data$aa],
           region = factor(region, c("WATL", "EATL")),
           Fleet = ifelse(ff == 1, "Otolith", "Genetic")) %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, Design$model_name)) %>%
  mutate(lwr = plogis(qlogis(obs) - 1.96 * se),
         upr = plogis(qlogis(obs) + 1.96 * se))

#g <- soo %>%
#  filter(stock == "WBFT") %>%
#  filter(grepl("SOO1", model)) %>%
#  filter(!is.na(pred), pred > 0) %>%
#  mutate(Age = paste("Age", Age)) %>%
#  ggplot(aes(year, obs, colour = Fleet, fill = Fleet)) +
#  geom_line(aes(y = pred, linetype = model), linewidth = 0.15, colour = "black") +
#  #geom_line(aes(y = pred), colour = 'red') +
#  geom_point(size = 0.5, shape = 1) +
#  geom_linerange(linewidth = 0.25, aes(ymin = lwr, ymax = upr)) +
#  facet_grid(vars(region), vars(Age)) +
#  coord_cartesian(ylim = c(0, 1)) +
#  #scale_shape_manual(values = c(1, 8)) +
#  #coord_cartesian(xlim = c(1970, 2025)) +
#  labs(x = "Year", y = "Proportion WBFT",
#       linetype = "Model",
#       fill = NULL, shape = NULL,
#       colour = NULL,
#       title = "Stock of origin (Set 1)") +
#  theme(legend.position = "bottom")
#ggsave("figures/fit/compare_SOO1_fit.png", g, height = 5, width = 8)
#
#g <- soo %>%
#  filter(stock == "WBFT") %>%
#  filter(grepl("SOO2", model)) %>%
#  filter(!is.na(pred), pred > 0) %>%
#  mutate(Age = paste("Age", Age)) %>%
#  ggplot(aes(year, obs, colour = Fleet, fill = Fleet)) +
#  geom_line(aes(y = pred, linetype = model), linewidth = 0.15, colour = "black") +
#  #geom_line(aes(y = pred), colour = 'red') +
#  geom_point(size = 0.5, shape = 1) +
#  geom_linerange(linewidth = 0.25, aes(ymin = lwr, ymax = upr)) +
#  facet_grid(vars(region), vars(Age)) +
#  coord_cartesian(ylim = c(0, 1)) +
#  #scale_shape_manual(values = c(1, 8)) +
#  #coord_cartesian(xlim = c(1970, 2025)) +
#  labs(x = "Year", y = "Proportion WBFT",
#       linetype = "Model",
#       fill = NULL, shape = NULL,
#       colour = NULL,
#       title = "Stock of origin (Set 2)") +
#  theme(legend.position = "bottom")
#ggsave("figures/fit/compare_SOO2_fit.png", g, height = 5, width = 8)

soo1_pred <- soo %>%
  filter(stock == "WBFT") %>%
  filter(model %in% Design$model_name[1:2]) %>%
  #filter(grepl("SOO1", model)) %>%
  mutate(Age = paste("Age", Age)) %>%
  filter(!is.na(pred), pred > 0)
soo1_obs <- soo %>%
  filter(stock == "WBFT") %>%
  filter(model %in% Design$model_name[1:2]) %>%
  #filter(grepl("SOO1", model)) %>%
  mutate(Age = paste("Age", Age)) %>%
  filter(!is.na(obs))

g <- soo1_pred %>%
  ggplot(aes(year, pred)) +
  geom_line(aes(colour = model), linewidth = 0.15) +
  geom_linerange(linewidth = 0.25, aes(ymin = lwr, ymax = upr)) +
  geom_point(data = soo1_obs, aes(y = obs, fill = Fleet), shape = 21, size = 0.75) +
  geom_line(data = soo1_obs, aes(y = obs, linetype = Fleet), linewidth = 0.1) +
  facet_grid(vars(region), vars(Age)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("black", "white")) +
  labs(x = "Year", y = "Proportion WBFT",
       colour = "Model prediction",
       linetype = "Data", fill = "Data",
       title = "Stock of origin (Set 1)") +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_SOO1_fit.png", g, height = 5, width = 8)


soo2_pred <- soo %>%
  filter(stock == "WBFT") %>%
  #filter(model %in% Design$model_name[3:4]) %>%
  filter(grepl("SOO2", model)) %>%
  mutate(Age = paste("Age", Age)) %>%
  filter(!is.na(pred), pred > 0)
soo2_obs <- soo %>%
  filter(stock == "WBFT") %>%
  #filter(model %in% Design$model_name[3:4]) %>%
  filter(grepl("SOO2", model)) %>%
  mutate(Age = paste("Age", Age)) %>%
  filter(!is.na(obs))

g <- soo2_pred %>%
  ggplot(aes(year, pred)) +
  geom_line(aes(colour = model), linewidth = 0.15) +
  geom_linerange(linewidth = 0.25, aes(ymin = lwr, ymax = upr)) +
  geom_point(data = soo2_obs, aes(y = obs, fill = Fleet), shape = 21, size = 0.75) +
  geom_line(data = soo2_obs, aes(y = obs, linetype = Fleet), linewidth = 0.15) +
  facet_grid(vars(region), vars(Age)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("black", "white")) +
  labs(x = "Year", y = "Proportion WBFT",
       colour = "Model prediction",
       linetype = "Data", fill = "Data",
       title = "Stock of origin (Set 2)") +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_SOO2_fit.png", g, height = 5, width = 8)

# Tag transitions
tags <- lapply(1:length(fits), function(i) {

  lapply(1:3, function(ac) {
    lapply(1:2, function(s) {
      x <- plot_tagmov(fits[[i]], s = s, aa = ac, figure = FALSE)
      if (is.null(x)) data.frame() else {
        x %>%
          mutate(N = sum(obs), .by = c(stock, aa, from, season)) %>%
          mutate(obs = obs/sum(obs), .by = c(stock, aa, from, season))
      }
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(obs = ifelse(is.na(obs), 0, obs)) %>%
  mutate(model = factor(model, Design$model_name))

plot_tags <- function(tags, title = NULL, type = c("departure", "arrival")) {
  type <- match.arg(type)

  tags_plot <- tags %>%
    mutate(m = strsplit(as.character(season), "Season ") %>% sapply(getElement, 2) %>% as.integer()) %>%
    mutate(season_arrive = paste("Season", ifelse(m == 4, 1, m + 1))) %>%
    mutate(from_num = match(from, c("GOM", "WATL", "EATL", "MED"))) %>%
    mutate(to_num = match(to, c("GOM", "WATL", "EATL", "MED"))) %>%
    mutate(from_label = paste("Origin:", from) %>% factor(paste("Origin:", c("GOM", "WATL", "EATL", "MED"))),
           to_label = paste("Destination:", to) %>% factor(paste("Destination:", c("GOM", "WATL", "EATL", "MED"))))

  if (type == "departure") {
    Nsamp <- tags_plot %>%
      filter(model %in% unique(model)[1]) %>%
      summarise(N_from = unique(N), .by = c(season, from, from_label, from_num))

    g <- tags_plot %>%
      ggplot(aes(to_num, pred)) +
      facet_grid(vars(from_label), vars(season)) +
      geom_text(data = Nsamp, aes(label = paste("N =", N_from)), x = Inf, y = Inf, hjust = "inward", vjust = "inward") +
      geom_line(aes(y = obs), colour = "black") +
      geom_point(aes(y = obs), shape = 1, colour = "black") +
      geom_line(aes(colour = model), linewidth = 1) +
      scale_x_continuous(breaks = 1:4, labels = c("GOM", "WATL", "EATL", "MED")) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Destination", y = "Proportion (departure from origin)", colour = "Model", title = title) +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 2))
  } else {

    g <- tags_plot %>%
      ggplot(aes(from_num, pred)) +
      facet_grid(vars(to_label), vars(season_arrive)) +
      geom_line(aes(y = obs), colour = "black") +
      geom_point(aes(y = obs), shape = 1, colour = "black") +
      geom_line(aes(colour = model), linewidth = 1) +
      scale_x_continuous(breaks = 1:4, labels = c("GOM", "WATL", "EATL", "MED")) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Origin", y = "Proportion (arrival)", colour = "Model", title = title) +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 2))
  }
  g
}



g <- filter(tags, stock == "EBFT", aa == 1, model %in% Design$model_name[1:4]) %>%
  plot_tags(title = paste("EBFT, Age", aa[1]))
ggsave("figures/fit/compare_tag_EBFT_a1.png", g, height = 6, width = 6)

g <- filter(tags, stock == "EBFT", aa == 2, model %in% Design$model_name[1:4]) %>%
  plot_tags(title = paste("EBFT, Age", aa[2]))
ggsave("figures/fit/compare_tag_EBFT_a2.png", g, height = 6, width = 6)

g <- filter(tags, stock == "EBFT", aa == 3, model %in% Design$model_name[1:4]) %>%
  plot_tags(title = paste("EBFT, Age", aa[3]))
ggsave("figures/fit/compare_tag_EBFT_a3.png", g, height = 6, width = 6)

g <- filter(tags, stock == "WBFT", model %in% Design$model_name[1:4]) %>%
  plot_tags(title = "WBFT")
ggsave("figures/fit/compare_tag_WBFT.png", g, height = 6, width = 6)

# CAL
dat <- get_MSAdata(fits[[1]])
CAL <- lapply(1:length(fits), function(i) {

  dat <- get_MSAdata(fits[[i]])

  lapply(1:dat@Dfishery@nf, function(f) {
    lapply(1:dat@Dmodel@nr, function(r) {
      x <- plot_CAL(fits[[i]], f = f, r = r, do_mean = FALSE, figure = FALSE)
      if (is.null(x)) x <- data.frame()
      return(x)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(model = Design$model_name[i])

}) %>%
  bind_rows()

CAL_agg <- CAL %>%
  mutate(  # for each time step, fleet, and region, re-do N
    obs2 = N * obs/sum(obs, na.rm = TRUE),
    pred2 = N * pred/sum(pred, na.rm = TRUE),
    .by = c(year, fleet, region, model)
  ) %>%
  summarise(
    obs = sum(obs2, na.rm = TRUE),
    pred = sum(pred2, na.rm = TRUE),
    .by = c(fleet, lmid, model)
  ) %>%
  mutate(
    obs = obs/sum(obs),
    pred = pred/sum(pred),
    .by = c(fleet, model)
  ) %>%
  mutate(fleet = factor(fleet, dat@Dlabel@fleet)) %>%
  mutate(model = factor(model, Design$model_name))

g <- CAL_agg %>%
  filter(model == Design$model_name[1]) %>%
  ggplot(aes(lmid, obs)) +
  geom_area(fill = "grey80", linewidth = 0.1, colour = "black") +
  geom_point() +
  geom_line(data = filter(CAL_agg, model %in% Design$model_name[1:4]),
            aes(y = pred, colour = model), linewidth = 1.25) +
  facet_wrap(vars(fleet), ncol = 3, scales = "free_y") +
  labs(x = "Length Bin", y = "Proportion", colour = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_CAL_agg_fit.png", g, height = 8, width = 6)

# Plot movement
mov <- lapply(1:nrow(Design), function(i) {
  lapply(1:2, function(s) {
    plot_mov(fits[[i]], s = s, figure = FALSE)
  }) %>%
    bind_rows() %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(Destination = ifelse(Destination == "Equilibrium", "Eq.", Destination)) %>%
  filter(model %in% Design$model_name[1:4])

plot_mov2 <- function(mov) {
  areas <- c("GOM", "WATL", "EATL", "MED", "", "Eq.")
  val <- seq(0, 1, 0.01)
  cols <- hcl.colors(length(val), palette = "Peach", rev = TRUE) %>% structure(names = val)

  g <- mov %>%
    mutate(x = match(Destination, areas), yy = match(Origin, areas)) %>%
    ggplot(aes(x, yy, fill = round(proportion, 2) %>% factor())) +
    geom_tile(colour = "black") +
    geom_text(size = 2.25, aes(label = round(proportion, 2))) +
    facet_grid(vars(model), vars(Season)) +
    scale_fill_manual(values = cols) +
    guides(fill = "none") +
    scale_x_continuous(breaks = 1:length(areas), labels = areas) +
    scale_y_continuous(breaks = 1:4, labels = areas[1:4]) +
    coord_transform(reverse = "y", expand = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Destination", y = "Origin", fill = "Proportion")
  g
}

g <- mov %>%
  filter(stock == "EBFT") %>%
  plot_mov2() +
  labs(title = "EBFT")
ggsave("figures/fit/compare_mov_EBFT.png", g, height = 6, width = 6)

g <- mov %>%
  filter(stock == "WBFT") %>%
  plot_mov2() +
  labs(title = "WBFT")
ggsave("figures/fit/compare_mov_WBFT.png", g, height = 6, width = 6)


# Catch residual
Cresid <- lapply(1:length(fits), function(i) {

  residuals(fits[[i]], vars = "Cobs_ymfr") %>%
    reshape2::melt() %>%
    mutate(model = Design$model_name[i])

}) %>%
  bind_rows() %>%
  mutate(model = factor(model, Design$model_name)) %>%
  filter(!is.na(value))

g <- Cresid %>%
  filter(abs(value) < 1) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "grey80", linewidth = 0.1, colour = "black") +
  facet_wrap(vars(model), ncol = 3) +
  labs(x = "Catch residual")
ggsave(file.path(dir_save, "Cresid.png"), g, height = 4, width = 6)


# Depletion
S_S0 <- lapply(1:nrow(Design), function(i) {
  dat <- get_MSAdata(fits[[i]])
  S_ys <- apply(fits[[i]]@report$S_yrs, c(1, 3), sum)
  S0 <- fits[[i]]@report$SB0_s
  dep <- structure(t(S_ys)/S0, dimnames = list(Stock = dat@Dlabel@stock, Year = dat@Dlabel@year))

  reshape2::melt(dep, value.name = "dep") %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows()

g <- S_S0 %>%
  filter(model %in% Design$model_name[1:4]) %>%
  ggplot(aes(Year, dep, linetype = Stock)) +
  geom_line() +
  facet_wrap(vars(model)) +
  coord_cartesian(ylim = c(0, 1.5)) +
  expand_limits(y = 0) +
  labs(y = expression(S/S[0]))
ggsave("figures/fit/compare_depletion.png", g, height = 4, width = 5)

# SSB by season
SB_season <- lapply(1:nrow(Design), function (i) {
  fit <- fits[[i]]
  dat <- get_MSAdata(fit)

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
  reshape2::melt(S_ymrs, value.name = "S") %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows()

for (i in 1:nrow(Design)) {
  g <- SB_season %>%
    filter(model == Design$model_name[i]) %>%
    #filter(stock == "EBFT") %>%
    ggplot(aes(year, S, fill = stock)) +
    geom_col(width = 1) +
    facet_grid(vars(region), vars(season)) +
    labs(x = "Year", y = "Mature biomass", fill = "Stock", title = Design$model_name[i]) +
    scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("figures/fit/SSB_area_season_model", i, ".png"), g, height = 5, width = 6)
}

# Total biomass by season
B_season <- lapply(1:nrow(Design), function (i) {
  dat <- get_MSAdata(fits[[i]])
  B <- plot_B(fits[[i]], figure = FALSE) %>%
    mutate(season = 4 * (year - floor(year)) + 1,
           Season = dat@Dlabel@season[season],
           model = Design$model_name[i])
  return(B)
}) %>%
  bind_rows()

for (i in 1:nrow(Design)) {
  g <- B_season %>%
    filter(model == Design$model_name[i]) %>%
    mutate(year = floor(year)) %>%
    #filter(stock == "EBFT") %>%
    ggplot(aes(year, B, fill = stock)) +
    geom_col(width = 1) +
    facet_grid(vars(region), vars(Season)) +
    labs(x = "Year", y = "Total biomass", fill = "Stock", title = Design$model_name[i]) +
    scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("figures/fit/SSB_area_season_model", i, ".png"), g, height = 5, width = 6)
}

# Calculate regional exploitation rate
u <- lapply(1:nrow(Design), function (i) {
  fit <- fits[[i]]
  dat <- get_MSAdata(fit)
  CB_ymrs <- apply(fit@report$CB_ymfrs, c(1, 2, 4, 5), sum)
  B_ymrs <- fit@report$B_ymrs
  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
  U_ymrs <- CB_ymrs/B_ymrs
  U_ymrs[CB_ymrs < 1e-8] <- 0
  U_yrs <- multiSA:::collapse_yearseason(U_ymrs) %>%
    structure(dimnames = list(Year = Year, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock))
  reshape2::melt(U_yrs, value.name = "Ex") %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(Season = 4 * (Year - floor(Year)) + 1)

u_eq <- lapply(1:nrow(Design), function (i) {
  fit <- fits[[i]]
  dat <- get_MSAdata(fit)

  CB_mrs <- apply(fit@report$initCB_mfrs, c(1, 3, 4), sum)
  N_mars <- sapply2(1:dat@Dmodel@ns, function(s) fit@report$initNPR_mars[, , , s] * fit@report$initReq_s[s])
  B_mrs <- sapply2(1:dat@Dmodel@nr, function(r) N_mars[, , r, ] * dat@Dstock@swt_ymas[1, , , ]) %>%
    apply(c(1, 4, 3), sum)

  U_mrs <- CB_mrs/B_mrs
  U_mrs[CB_mrs < 1e-8] <- 0

  structure(U_mrs, dimnames = list(Season = dat@Dlabel@season, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock)) %>%
    reshape2::melt(value.name = "Ex") %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(Year = min(u$Year) - 1)

for (i in 1:nrow(Design)) {
  g <- u %>%
    filter(model == Design$model_name[i]) %>%
    mutate(Season = paste("Season", Season)) %>%
    ggplot(aes(floor(Year), Ex, colour = Stock)) +
    facet_grid(vars(Region), vars(Season), scales = "free") +
    geom_line() +
    geom_point(data = u_eq %>% filter(model == Design$model_name[i])) +
    #geom_point(alpha = 0.5, size = 0.75, aes(colour = factor(Season))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Year", y = "Seasonal Catch/Biomass", title = Design$model_name[i]) +
    theme(legend.position = "bottom")
  ggsave(paste0("figures/fit/regional_exploitation_model", i, ".png"), g, height = 5, width = 6)
}

