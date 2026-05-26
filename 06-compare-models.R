
library(multiSA)
library(tidyverse)

Design <- readr::read_csv("tables/Design_05.22.2026.csv")

fits <- lapply(Design$output_name, function(i) readRDS(file.path("model_output", paste0("fit_", i, ".rds"))))


# Compare SSB
SSB <- lapply(1:nrow(Design), function(i) {
  plot_S(fits[[i]], s = 1:2, figure = FALSE) %>%
    mutate(model = Design$model_name[i])
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, Design$model_name))

prior <- data.frame(
  model = Design$model_name[2],
  stock = "WBFT",
  year = 2021,
  S = 18000,
  lwr = exp(log(18000) - 1.96 * 0.18),
  upr = exp(log(18000) + 1.96 * 0.18)
) %>%
  mutate(model = factor(model, Design$model_name))

g <- SSB %>%
  ggplot(aes(year, S)) +
  facet_grid(vars(stock), vars(model)) +
  geom_col(width = 1, aes(fill = region)) +
  geom_pointrange(data = prior, size = 0.25, aes(ymin = lwr, ymax = upr)) +
  labs(x = "Year", y = "Spawning stock biomass", fill = NULL) +
  scale_fill_manual(values = multiSA:::make_color(4, "region"))
ggsave("figures/fit/compare_SSB.png", g, height = 5, width = 7)


# Likelihoods
like <- sapply(fits, function(i) {
  sapply(i@report[grepl("loglike", names(i@report))], sum)
}) %>%
  `colnames<-`(Design$model_name)

like[, 1:2] %>%
  apply(1, function(x) x - max(x)) %>%
  t()


dat <- get_MSAdata(fit)

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
  theme(legend.position = "bottom")
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

g <- soo %>%
  filter(stock == "WBFT") %>%
  filter(!is.na(pred), pred > 0) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(year, obs, colour = Fleet, fill = Fleet)) +
  geom_line(aes(y = pred, linetype = model), linewidth = 0.15, colour = "black") +
  #geom_line(aes(y = pred), colour = 'red') +
  geom_point(size = 0.5, shape = 1) +
  geom_linerange(linewidth = 0.25, aes(ymin = lwr, ymax = upr)) +
  facet_grid(vars(region), vars(Age)) +
  coord_cartesian(ylim = c(0, 1)) +
  #scale_shape_manual(values = c(1, 8)) +
  #coord_cartesian(xlim = c(1970, 2025)) +
  labs(x = "Year", y = "Proportion WBFT",
       linetype = "Model",
       fill = NULL, shape = NULL,
       colour = NULL,
       title = "Stock of origin") +
  theme(legend.position = "bottom")
ggsave("figures/fit/compare_SOO_fit.png", g, height = 5, width = 8)

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
      summarise(N_from = unique(N), .by = c(season, from))

    g <- tags_plot %>%
      ggplot(aes(to_num, pred)) +
      facet_grid(vars(from_label), vars(season)) +
      geom_label(data = Nsamp, aes(label = paste("N =", N_from)), x = Inf, y = Inf, hjust = "inward", vjust = "inward") +
      geom_line(aes(y = obs), colour = "black") +
      geom_point(aes(y = obs), shape = 1, colour = "black") +
      geom_line(aes(colour = model), linewidth = 1) +
      scale_x_continuous(breaks = 1:4, labels = c("GOM", "WATL", "EATL", "MED")) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Destination", y = "Proportion (departure)", colour = "Model", title = title) +
      theme(legend.position = "bottom")
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
      theme(legend.position = "bottom")
  }
  g
}

g <- filter(tags, stock == "WBFT") %>%
  plot_tags(title = "WBFT")

g <- filter(tags, stock == "WBFT") %>%
  plot_tags(title = "WBFT", type = "arrival")

g <- filter(tags, stock == "EBFT", aa == 1) %>%
  plot_tags(title = paste("EBFT, Age", aa[1]))

g <- filter(tags, stock == "EBFT", aa == 2) %>%
  plot_tags(title = paste("EBFT, Age", aa[2]))

g <- filter(tags, stock == "EBFT", aa == 3) %>%
  plot_tags(title = paste("EBFT, Age", aa[3]))


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
  geom_line(data = CAL_agg, aes(y = pred, colour = model), linewidth = 1.25) +
  facet_wrap(vars(fleet), ncol = 3, scales = "free_y") +
  labs(x = "Length Bin", y = "Proportion", colour = NULL)
ggsave(file.path(dir_save, "CAL_agg_fit.png"), g, height = 8, width = 6)

