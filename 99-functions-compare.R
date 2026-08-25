
library(tidyverse)
library(multiSA)

# Diagnostics with hessian matrix
get_corr <- function(fit) {
  corr <- fit@SD$cov.fixed %>%
    structure(dimnames = dimnames(fit@SD$env$hessian)) %>%
    cov2cor() %>%
    round(4)

  corr[upper.tri(corr)] <- NA

  reshape2::melt(corr) %>%
    filter(Var1 != Var2) %>%
    arrange(desc(abs(value)))
}

condition_number <- function(fit) {
  ev <- eigen(fit@SD$env$hessian)
  ev$values[1]/ev$values[length(ev$values)]
}

# Francis weights
Francis_weights <- function(fit) {
  dat <- get_MSAdata(fit)

  # Annual model only
  mlen <- lapply(1:dat@Dfishery@nf, function(f) {
    lapply(1:dat@Dmodel@nr, function(r) {
      x <- plot_CAL(fit, f = f, r = r, do_mean = TRUE, figure = FALSE)
      if (is.null(x)) x <- data.frame()
      return(x)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    rename(mpred = pred, mobs = obs)

  CAL <- lapply(1:dat@Dfishery@nf, function(f) {
    lapply(1:dat@Dmodel@nr, function(r) {
      x <- plot_CAL(fit, f = f, r = r, do_mean = FALSE, figure = FALSE)
      if (is.null(x)) x <- data.frame()
      return(x)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()

  wt <- CAL %>%
    filter(N > 0) %>%
    mutate(Npred = pred * N) %>%
    left_join(mlen, by = c("year", "fleet", "region")) %>%
    mutate(denom = sum(Npred), .by = c(year, fleet)) %>%
    mutate(phat_y = Npred/denom) %>%
    summarise(VAR_num = sum(phat_y * (lmid - mpred)^2), N = unique(N), .by = c(year, fleet, region)) %>%
    left_join(mlen, by = c("year", "fleet", "region")) %>%
    mutate(SE = sqrt(VAR_num/N), Z = (mobs - mpred)/SE) %>%
    summarise(w = 1/var(Z, na.rm = TRUE), .by = fleet)
  return(wt)
}


# Plotting functions
plot_SSB <- function(fits, model_name, scales = "fixed") {
  SSB <- lapply(1:length(model_name), function(i) {
    plot_S(fits[[i]], s = 1:2, figure = FALSE) %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name))

  prior <- data.frame(
    model = model_name,
    stock = "WBFT",
    year = 2018,
    S = 22000,
    lwr = exp(log(22000) - 1.96 * 0.18),
    upr = exp(log(22000) + 1.96 * 0.18)
  ) %>%
    mutate(model = factor(model, model_name))

  g <- SSB %>%
    ggplot(aes(year, S)) +
    facet_grid(vars(stock), vars(model), scales = scales) +
    geom_col(width = 1, aes(fill = region)) +
    geom_pointrange(data = prior, size = 0.25, aes(ymin = lwr, ymax = upr)) +
    labs(x = "Year", y = "Spawning stock biomass (season 2)", fill = NULL) +
    scale_fill_manual(values = multiSA:::make_color(length(unique(SSB$region)), "region")) +
    theme(legend.position = "bottom")

  g
}

plot_rec <- function(fits, model_name) {
  rec <- lapply(1:length(model_name), function(i) {
    plot_R(fits[[i]], figure = FALSE) %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name))

  g <- rec %>%
    mutate(year = as.numeric(year)) %>%
    ggplot(aes(year, R)) +
    facet_grid(vars(stock), vars(model), scales = "free_y") +
    geom_line() +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Recruitment")

  g
}

plot_recdev <- function(fits, model_name) {
  recdev <- lapply(1:length(model_name), function(i) {
    lapply(1:2, function(s) {
      plot_Rdev(fits[[i]], figure = FALSE, s = s) %>%
        mutate(model = model_name[i])
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name))

  g <- recdev %>%
    mutate(year = as.numeric(year)) %>%
    ggplot(aes(year, dev)) +
    geom_line(linewidth = 0.25) +
    geom_point(size = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(vars(stock), vars(model), scales = "free_y") +
    labs(x = "Year", y = "Recruitment deviation")

  g
}

plot_sel <- function(fits, model_name) {

  dat <- get_MSAdata(fits[[1]])
  sel <- lapply(1:length(model_name), function(i) {
    dat <- get_MSAdata(fits[[i]])
    lapply(1:dat@Dfishery@nf, function(f) {
      plot_self(fits[[i]], f = f, figure = FALSE)
    }) %>%
      bind_rows() %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name),
           fleet = factor(fleet, dat@Dlabel@fleet))

  g <- sel %>%
    ggplot(aes(length, sel, colour = model, linetype = year)) +
    facet_wrap(vars(fleet), ncol = 3) +
    geom_line() +
    labs(x = "Length", y = "Selectivity", linetype = "Selectivity\nblock", colour = "Model") +
    scale_linetype_manual(values = c(1, 3, 2)) +
    guides(linetype = guide_legend(ncol = 1), colour = guide_legend(ncol = 1)) +
    theme(legend.position = "bottom")

  g
}

plot_sel_index <- function(fits, model_name, type = c("age", "length")) {
  type <- match.arg(type)

  dat <- get_MSAdata(fits[[1]])
  seli <- lapply(1:length(model_name), function(i) {
    lapply(1:dat@Dsurvey@ni, function(ii) {
      plot_seli(fits[[i]], i = ii, figure = FALSE)
    }) %>%
      bind_rows() %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name),
           name = factor(name, dat@Dlabel@index))
  if (type == "length") {
    g <- seli %>%
      filter(!is.na(length)) %>%
      ggplot(aes(length, sel, colour = model)) +
      facet_wrap(vars(name), ncol = 3) +
      geom_line() +
      labs(x = "Length", y = "Selectivity", colour = NULL) +
      theme(legend.position = "bottom")
  } else {
    g <- seli %>%
      filter(is.na(length)) %>%
      ggplot(aes(age, sel)) +
      facet_wrap(vars(name), ncol = 3) +
      geom_line() +
      labs(x = "Age", y = "Selectivity") +
      theme(legend.position = "bottom")
  }

  g
}

plot_fit_index <- function(fits, model_name, type = c("cpue", "fi")) {
  type <- match.arg(type)
  dat <- get_MSAdata(fits[[1]])

  index_all <- lapply(1:length(fits), function(i) {
    plot_index(fits[[i]], i = 1:dat@Dsurvey@ni, zoom = TRUE, figure = FALSE) %>%
      mutate(resid = log(obs/pred)) %>%
      mutate(name = factor(name, dat@Dlabel@index)) %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, Design$model_name)) %>%
    filter(!is.na(name))

  if (type == "cpue") {

    g <- index_all %>%
      filter(name %in% dat@Dlabel@index[1:16]) %>%
      ggplot(aes(year, obs)) +
      geom_point(size = 0.5) +
      geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr)) +
      geom_line(aes(y = pred, colour = model)) +
      expand_limits(y = 0) +
      facet_wrap(vars(name), ncol = 3, scales = "free_y") +
      labs(x = "Year", y = "CPUE", colour = "Model") +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 2))
  } else {
    g <- index_all %>%
      filter(!name %in% dat@Dlabel@index[1:16]) %>%
      ggplot(aes(year, obs)) +
      geom_point(size = 0.5) +
      geom_linerange(linewidth = 0.1, aes(ymin = lwr, ymax = upr)) +
      geom_line(aes(y = pred, colour = model)) +
      expand_limits(y = 0) +
      facet_wrap(vars(name), ncol = 3, scales = "free_y") +
      labs(x = "Year", y = "Index", colour = "Model") +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 2))
  }

  g
}


plot_fit_soo3 <- function(fits, model_name) {
  soo3 <- lapply(1:length(fits), function(i) {
    dat <- get_MSAdata(fits[[i]])

    lapply(1:nrow(dat@Dfishery@SC_ff), function(ff) {
      fleet_name <- dat@Dlabel@fleet[as.logical(dat@Dfishery@SC_ff[ff, ])]
      plot_SC(fits[[i]], r = ifelse(dat@Dmodel@nr == 4, 2, 1), ff = ff, prop = TRUE, figure = FALSE) %>%
        mutate(fleet = fleet_name)
    }) %>%
      bind_rows() %>%
      mutate(region = "WATL") %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(model = factor(model, model_name)) %>%
    mutate(lwr = plogis(qlogis(obs) - 1.96 * se),
           upr = plogis(qlogis(obs) + 1.96 * se)) %>%
    mutate(N = sum(pred, na.rm = TRUE), .by = c(year, ff, aa, region, model)) %>%
    filter(N > 0)

  g <- soo3 %>%
    mutate(season = 4 * (year - floor(year)) + 1) %>%
    filter(stock == "WBFT", year < 2025, season > 2) %>%
    ggplot(aes(year, pred)) +
    geom_line(aes(colour = model)) +
    geom_pointrange(linewidth = 0.25, size = 0.1, aes(y = obs, ymin = lwr, ymax = upr)) +
    facet_grid(vars(fleet), vars(model)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Year", y = "Proportion WBFT",
         colour = "Model prediction",
         title = "Stock mixing in WATL") +
    theme(legend.position = "bottom")

  g
}

plot_fit_CAL_agg <- function(fits, model_name) {

  dat <- get_MSAdata(fits[[1]])

  CAL <- lapply(1:length(fits), function(i) {
    lapply(1:dat@Dfishery@nf, function(f) {
      lapply(1:dat@Dmodel@nr, function(r) {
        x <- plot_CAL(fits[[i]], f = f, r = r, do_mean = FALSE, agg = FALSE, figure = FALSE)
        if (is.null(x)) x <- data.frame()
        return(x)
      }) %>%
        bind_rows()
    }) %>%
      bind_rows() %>%
      mutate(model = model_name[i])

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
    mutate(model = factor(model, model_name))

  g <- CAL_agg %>%
    ggplot(aes(lmid)) +
    geom_area(data = filter(CAL_agg, model == model_name[1]),
              aes(y = obs), fill = "grey80", linewidth = 0.1, colour = "black") +
    geom_point(data = filter(CAL_agg, model == model_name[1]), aes(y = obs)) +
    geom_line(aes(y = pred, colour = model), linewidth = 1.25) +
    facet_wrap(vars(fleet), ncol = 3, scales = "free_y") +
    labs(x = "Length Bin", y = "Proportion", colour = NULL) +
    theme(legend.position = "bottom")

  g
}

plot_mlen <- function(fits, model_name) {

  dat <- get_MSAdata(fits[[1]])

  panel_factor <- outer(dat@Dlabel@fleet, dat@Dlabel@region, function(i, j) paste0(i, " (", j, ")")) %>%
    t() %>%
    as.character()

  mlen <- lapply(1:length(fits), function(i) {
    plot_CAL(fits[[i]], do_mean = TRUE, figure = FALSE) %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(fleet = factor(fleet, dat@Dlabel@fleet)) %>%
    mutate(region = factor(region, dat@Dlabel@region)) %>%
    arrange(year) %>%
    mutate(panel = paste0(fleet, " (", region, ")") %>% factor(panel_factor))

  mlen_obs <- filter(mlen, !is.na(obs))
  mlen_pred <- filter(mlen, !is.na(pred))
  g <- mlen_obs %>%
    ggplot(aes(year, obs)) +
    geom_point(shape = 1, size = 0.75) +
    geom_line(linewidth = 0.1) +
    geom_line(data = mlen_pred, aes(y = pred, colour = model)) +
    facet_wrap(vars(panel), ncol = 4) +
    theme(legend.position = "bottom") +
    labs(x = "Year", y = "Mean length", fill = NULL, colour = NULL)

  g
}



## Tag transitions ----
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

plot_tags <- function(fits, model_name, title = NULL, type = c("departure", "arrival"),
                      stock = c("EBFT", "WBFT"), ac = 1) {
  stock <- match.arg(stock)
  type <- match.arg(type)

  tags <- lapply(1:length(fits), function(i) {
    x <- plot_tagmov(fits[[i]], s = ifelse(stock == "EBFT", 1, 2), aa = ac, figure = FALSE)
    if (is.null(x)) {
      data.frame()
    } else {
      x %>%
        mutate(N = sum(obs), .by = c(stock, aa, from, season)) %>%
        mutate(obs = obs/sum(obs), .by = c(stock, aa, from, season)) %>%
        mutate(model = model_name[i])
    }
  }) %>%
    bind_rows() %>%
    mutate(obs = ifelse(is.na(obs), 0, obs)) %>%
    mutate(model = factor(model, model_name))

  tags_plot <- tags %>%
    mutate(m = strsplit(as.character(season), "Season ") %>% sapply(getElement, 2) %>% as.integer()) %>%
    mutate(season_arrive = paste("Season", ifelse(m == 4, 1, m + 1))) %>%
    mutate(from_num = match(from, c("GOM", "WATL", "EATL", "MED"))) %>%
    mutate(to_num = match(to, c("GOM", "WATL", "EATL", "MED"))) %>%
    mutate(from_label = paste("Origin:", from) %>% factor(paste("Origin:", c("GOM", "WATL", "EATL", "MED"))),
           to_label = paste("Destination:", to) %>% factor(paste("Destination:", c("GOM", "WATL", "EATL", "MED"))))

  if (type == "departure") {
    Nsamp <- tags_plot %>%
      filter(model == model_name[1]) %>%
      mutate(N_from = unique(N), .by = c(season, from, from_label, from_num))

    g <- tags_plot %>%
      ggplot(aes(to_num, pred)) +
      facet_grid(vars(from_label), vars(season)) +
      geom_text(data = Nsamp, aes(label = paste("N =", N_from)),
                x = ifelse(stock == "EBFT", -Inf, Inf), y = 1, hjust = "inward", vjust = "inward") +
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


plot_mov_gg <- function(fits, model_name, stock = c("EBFT", "WBFT")) {
  stock <- match.arg(stock)

  mov <- lapply(1:length(fits), function(i) {
    plot_mov(fits[[i]], s = ifelse(stock == "EBFT", 1, 2), figure = FALSE) %>%
      mutate(model = model_name[i])
  }) %>%
    bind_rows() %>%
    mutate(Destination = ifelse(Destination == "Equilibrium", "Eq.", Destination))

  areas <- c("GOM", "WATL", "EATL", "MED", "", "Eq.")
  val <- seq(0, 1, 0.01)
  cols <- hcl.colors(length(val), palette = "Peach", rev = TRUE) %>% structure(names = val)

  unique_areas <- unique(mov$Origin)
  areas <- c(unique_areas, "", "Eq.")

  g <- mov %>%
    mutate(x = match(Destination, areas), yy = match(Origin, areas)) %>%
    ggplot(aes(x, yy, fill = round(proportion, 2) %>% factor())) +
    geom_tile(colour = "black") +
    geom_text(size = 2.75, aes(label = round(proportion, 2) %>% format())) +
    facet_grid(vars(model), vars(Season)) +
    scale_fill_manual(values = cols) +
    guides(fill = "none") +
    scale_x_continuous(breaks = 1:length(areas), labels = areas) +
    scale_y_continuous(breaks = 1:length(unique_areas), labels = areas[1:length(unique_areas)]) +
    coord_transform(reverse = "y", expand = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Destination", y = "Origin", fill = "Proportion")
  g
}

## SRR ----
#SRR <- lapply(1:nrow(Design), function(i) {
#  lapply(1:2, function(s) {
#    fit <- fits[[i]]
#    dat <- get_MSAdata(fit)
#    Dlabel <- get_MSAdata(fit)@Dlabel
#
#    S_y <- apply(fit@report$S_yrs[, , s, drop = FALSE], 1, sum)
#    R_y <- fit@report$R_ys[, s]
#
#    Rpred_y <- R_y/fit@report$Rdev_ys[, s]
#
#    a <- fit@report$sralpha_s[s]
#    b <- fit@report$srbeta_s[s]
#
#    S_SRR <- seq(0, max(S_y), length.out = 100)
#    R_SRR <- calc_recruitment(S_SRR, SRR = dat@Dstock@SRR_s[s], a = a, b = b)
#
#    S2 <- S_y[-1]
#    R2 <- Rpred_y[-1]
#
#    list(
#      hist = data.frame(year = dat@Dlabel@year, S = S_y, R = R_y,
#                        model = Design$model_name[i],
#                        stock = dat@Dlabel@stock[s]),
#      pred = data.frame(S = S_SRR, R = R_SRR, model = Design$model_name[i],
#                        stock = dat@Dlabel@stock[s]),
#      phi = data.frame(phi = fit@report$phi_s[s], model = Design$model_name[i],
#                       stock = dat@Dlabel@stock[s])
#    )
#  })
#})
#
#SRR_hist <- lapply(SRR, function(i) lapply(i, getElement, 'hist') %>%
#                     bind_rows()) %>%
#  bind_rows()
#SRR_pred <- lapply(SRR, function(i) lapply(i, getElement, 'pred') %>%
#                     bind_rows()) %>%
#  bind_rows()
#
#SRR_phi <- lapply(SRR, function(i) lapply(i, getElement, 'phi') %>%
#                    bind_rows()) %>%
#  bind_rows()
#
#g <- ggplot(SRR_hist, aes(S, R)) +
#  geom_line(data = SRR_pred) +
#  geom_abline(data = SRR_phi, aes(intercept = 0, slope = 1/phi), linetype = 2) +
#  geom_point(size = 1, stroke = 0.1, shape = 21, aes(fill = year)) +
#  facet_grid(vars(stock), vars(model)) +
#  scale_fill_viridis_c() +
#  theme(legend.position = "bottom") +
#  labs(x = "Spawning biomass", y = "Recruitment", fill = "Year")
#ggsave(file.path(dir_save, "compare_SRR.png"), g, height = 4, width = 6)
#
#
#
#
## Catch residual
#Cresid <- lapply(1:length(fits), function(i) {
#
#  residuals(fits[[i]], vars = "Cobs_ymfr") %>%
#    reshape2::melt() %>%
#    mutate(model = Design$model_name[i])
#
#}) %>%
#  bind_rows() %>%
#  mutate(model = factor(model, Design$model_name)) %>%
#  filter(!is.na(value))
#
#g <- Cresid %>%
#  filter(abs(value) < 1) %>%
#  ggplot(aes(x = value)) +
#  geom_histogram(fill = "grey80", linewidth = 0.1, colour = "black") +
#  facet_wrap(vars(model), ncol = 3) +
#  labs(x = "Catch residual")
#
## Depletion ----
#S_S0 <- lapply(1:nrow(Design), function(i) {
#  dat <- get_MSAdata(fits[[i]])
#  S_ys <- apply(fits[[i]]@report$S_yrs, c(1, 3), sum)
#  S0 <- fits[[i]]@report$SB0_s
#  dep <- structure(t(S_ys)/S0, dimnames = list(Stock = dat@Dlabel@stock, Year = dat@Dlabel@year))
#
#  reshape2::melt(dep, value.name = "dep") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows()
#
#g <- S_S0 %>%
#  filter(model %in% Design$model_name[1:4]) %>%
#  ggplot(aes(Year, dep, linetype = Stock)) +
#  geom_line() +
#  facet_wrap(vars(model)) +
#  coord_cartesian(ylim = c(0, 1.5)) +
#  expand_limits(y = 0) +
#  labs(y = expression(S/S[0]))
#ggsave(file.path(dir_save, "compare_depletion.png"), g, height = 3.5, width = 5)
#
## SSB by season ----
plot_SB_season <- function(fit) {
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
  SB_season <- reshape2::melt(S_ymrs, value.name = "S")

  g <- SB_season %>%
    #filter(stock == "EBFT") %>%
    ggplot(aes(year, S, fill = stock)) +
    geom_col(width = 1) +
    facet_grid(vars(region), vars(season)) +
    labs(x = "Year", y = "Mature biomass", fill = "Stock") +
    scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  g
}


## Total biomass by season ----
#if (FALSE) {
#  B_season <- lapply(1:nrow(Design), function (i) {
#    dat <- get_MSAdata(fits[[i]])
#    B <- plot_B(fits[[i]], figure = FALSE) %>%
#      mutate(season = 4 * (year - floor(year)) + 1,
#             Season = dat@Dlabel@season[season],
#             model = Design$model_name[i])
#    return(B)
#  }) %>%
#    bind_rows()
#
#  for (i in 1:nrow(Design)) {
#    g <- B_season %>%
#      filter(model == Design$model_name[i]) %>%
#      mutate(year = floor(year)) %>%
#      #filter(stock == "EBFT") %>%
#      ggplot(aes(year, B, fill = stock)) +
#      geom_col(width = 1) +
#      facet_grid(vars(region), vars(Season)) +
#      labs(x = "Year", y = "Total biomass", fill = "Stock", title = Design$model_name[i]) +
#      scale_fill_manual(values = grDevices::hcl.colors(2, palette = "Set2")) +
#      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
#    ggsave(file.path(dir_save, paste0("B_area_season_model", i, ".png")), g, height = 5, width = 6)
#  }
#}
#
## Calculate seasonal, spatial exploitation rates ----
#u <- lapply(1:nrow(Design), function (i) {
#  fit <- fits[[i]]
#  dat <- get_MSAdata(fit)
#  CB_ymrs <- apply(fit@report$CB_ymfrs, c(1, 2, 4, 5), sum)
#  B_ymrs <- fit@report$B_ymrs
#  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
#  U_ymrs <- CB_ymrs/B_ymrs
#  U_ymrs[CB_ymrs < 1e-8] <- 0
#  U_yrs <- multiSA:::collapse_yearseason(U_ymrs) %>%
#    structure(dimnames = list(Year = Year, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock))
#  reshape2::melt(U_yrs, value.name = "Ex") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows() %>%
#  mutate(Season = 4 * (Year - floor(Year)) + 1)
#
#u_eq <- lapply(1:nrow(Design), function (i) {
#  fit <- fits[[i]]
#  dat <- get_MSAdata(fit)
#
#  CB_mrs <- apply(fit@report$initCB_mfrs, c(1, 3, 4), sum)
#  N_mars <- sapply2(1:dat@Dmodel@ns, function(s) fit@report$initNPR_mars[, , , s] * fit@report$initReq_s[s])
#  B_mrs <- sapply2(1:dat@Dmodel@nr, function(r) N_mars[, , r, ] * dat@Dstock@swt_ymas[1, , , ]) %>%
#    apply(c(1, 4, 3), sum)
#
#  U_mrs <- CB_mrs/B_mrs
#  U_mrs[CB_mrs < 1e-8] <- 0
#
#  structure(U_mrs, dimnames = list(Season = dat@Dlabel@season, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock)) %>%
#    reshape2::melt(value.name = "Ex") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows() %>%
#  mutate(Year = min(u$Year) - 1)
#
#for (i in 1:nrow(Design)) {
#  g <- u %>%
#    filter(model == Design$model_name[i]) %>%
#    mutate(Season = paste("Season", Season)) %>%
#    ggplot(aes(floor(Year), Ex, colour = Stock)) +
#    facet_grid(vars(Region), vars(Season), scales = "free") +
#    geom_line() +
#    geom_point(data = u_eq %>% filter(model == Design$model_name[i])) +
#    #geom_point(alpha = 0.5, size = 0.75, aes(colour = factor(Season))) +
#    coord_cartesian(ylim = c(0, 1)) +
#    labs(x = "Year", y = "Seasonal Catch/Biomass", title = Design$model_name[i]) +
#    theme(legend.position = "bottom")
#  ggsave(file.path(dir_save, paste0("spatial_exploitation_model", i, ".png")), g, height = 5, width = 6)
#}
#
## Calculate seasonal, spatial F ----
#FF <- lapply(1:nrow(Design), function (i) {
#  fit <- fits[[i]]
#  dat <- get_MSAdata(fit)
#  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
#
#  F_yrs <- apply(fit@report$F_ymars, c(1, 2, 4, 5), max) %>%
#    multiSA:::collapse_yearseason() %>%
#    structure(dimnames = list(Year = Year, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock))
#  reshape2::melt(F_yrs, value.name = "FM") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows() %>%
#  mutate(Season = 4 * (Year - floor(Year)) + 1)
#
#for (i in 1:nrow(Design)) {
#  g <- FF %>%
#    filter(model == Design$model_name[i]) %>%
#    mutate(Season = paste("Season", Season)) %>%
#    ggplot(aes(floor(Year), FM, colour = Stock)) +
#    facet_grid(vars(Region), vars(Season), scales = "free") +
#    geom_line() +
#    geom_jitter() +
#    expand_limits(y = 0) +
#    #coord_cartesian(ylim = c(0, 1)) +
#    labs(x = "Year", y = "Seasonal Apical F", title = Design$model_name[i]) +
#    theme(legend.position = "bottom")
#  ggsave(file.path(dir_save, paste0("spatial_F_model", i, ".png")), g, height = 5, width = 6)
#}
#
#
## Calculate seasonal, spatial F at age ----
## F is identical between stocks, what matters is availability (whether the fish are there or not)
#FF_age <- lapply(1:nrow(Design), function (i) {
#  fit <- fits[[i]]
#  dat <- get_MSAdata(fit)
#  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
#
#  F_yars <- multiSA:::collapse_yearseason(fit@report$F_ymars) %>%
#    structure(dimnames = list(Year = Year, Age = dat@Dlabel@age, Region = dat@Dlabel@region, Stock = dat@Dlabel@stock))
#
#  reshape2::melt(F_yars, value.name = "FM") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows() %>%
#  mutate(Season = 4 * (Year - floor(Year)) + 1)
#
#for (i in 1:nrow(Design)) {
#  g <- FF_age %>%
#    filter(Stock == "WBFT") %>%
#    filter(model == Design$model_name[i]) %>%
#    mutate(Season = paste("Season", Season)) %>%
#    ggplot(aes(floor(Year), Age, fill = FM)) +
#    facet_grid(vars(Region), vars(Season), scales = "free") +
#    geom_tile() +
#    scale_fill_viridis_c(option = "C") +
#    labs(x = "Year", y = "Age", fill = "Fishing mortality", title = Design$model_name[i]) +
#    theme(legend.position = "bottom")
#  ggsave(file.path(dir_save, paste0("spatial_F_model", i, ".png")), g, height = 5, width = 6)
#}
#
## Catch at age - aggregate across both stocks ----
#CAA <- lapply(1:nrow(Design), function (i) {
#  fit <- fits[[i]]
#  dat <- get_MSAdata(fit)
#  Year <- multiSA:::make_yearseason(dat@Dlabel@year, 4)
#
#  C_yar <- fit@report$CN_ymafrs %>%
#    apply(c(1, 2, 3, 5), sum) %>%
#    multiSA:::collapse_yearseason() %>%
#    structure(dimnames = list(Year = Year, Age = dat@Dlabel@age, Region = dat@Dlabel@region))
#
#  reshape2::melt(C_yar, value.name = "value") %>%
#    mutate(model = Design$model_name[i])
#}) %>%
#  bind_rows() %>%
#  mutate(Season = 4 * (Year - floor(Year)) + 1)
#
#for (i in 1:nrow(Design)) {
#  g <- CAA %>%
#    filter(model == Design$model_name[i]) %>%
#    mutate(Season = paste("Season", Season)) %>%
#    ggplot(aes(floor(Year), Age, fill = value)) +
#    facet_grid(vars(Region), vars(Season), scales = "free") +
#    geom_tile() +
#    scale_fill_viridis_c(option = "C") +
#    labs(x = "Year", y = "Age", fill = "Catch at age", title = Design$model_name[i]) +
#    theme(legend.position = "bottom")
#  ggsave(file.path(dir_save, paste0("spatial_CAA_model", i, ".png")), g, height = 5, width = 6)
#}


