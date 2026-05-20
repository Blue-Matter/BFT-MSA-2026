
library(multiSA)
library(tidyverse)

#### Master file to organize most model settings ----
xlsx_file <- file.path("data", "1_M3_data", "ICCAT_MSA_Data_2026Apr_v1.xlsx")

#### Area names ----
area_names <- readxl::read_excel(xlsx_file, sheet = "Areas") %>%
  select(Area, Name) %>%
  filter(1:nrow(.) %in% seq(1, nrow(.), 2))

#### Age classes ----
ageclass_key <- readxl::read_excel(xlsx_file, sheet = "Age_classes") %>%
  mutate(Name = ifelse(Class == 1, "0-4", ifelse(Class == 2, "5-8", "9+"))) %>%
  summarise(n = n(), .by = c(Class, Name))

#### Fleet names ----
fleet_names <- readxl::read_excel(xlsx_file, sheet = "Fleets") %>%
  mutate(FleetName = paste0(Number, ": ", Code))

#### Model Settings (Years, Age, Areas, Stocks, Seasons, Length bins, Priors) ----
MetaData <- readxl::read_excel(xlsx_file, sheet = "Meta_Data")
ModelYear <- seq(1950, 2025)

ny <- length(ModelYear) # 76
nr <- 4  # areas
ns <- 2  # stocks
na <- 36 # ages: 0, 1, 2, ... 35
nm <- 4  # seasons

len_bin <- readxl::read_excel(xlsx_file, sheet = "Length_classes")
nl <- nrow(len_bin) # 16

# MSA object for model structure
Dmodel <- new(
  "Dmodel",
  ny = ny,
  nm = nm,
  na = na,
  nl = nl,
  nr = nr,
  ns = ns,
  lbin = len_bin$LengthClass, # Length-16
  lmid = len_bin$LengthClass + 0.5 * unique(diff(len_bin$LengthClass)),
  Fmax = 3,
  y_phi = 1,
  scale_s = c(1, 10), # Multiply WBFT R0 by 10 to get to similar order of magnitude as EBFT to aid estimation
  nyinit = 2 * na, # Spool-up for 2 life cycles
  condition = "catch"
)

# Add priors here as a character string that is evaluated during the model run
# It really helps to know the internal model code in multiSA:::.MSA()

# Constrain rec devs as devvector (sums to zero)
prior_recdev <- paste0("dnorm(sum(p$log_rdev_ys[, ", 1:2, "]), 0, 0.01, log = TRUE)")

# Prior for proportion of stock in natal region in each season
prior_dist <- local({
  d <- readxl::read_excel(xlsx_file, sheet = "Seasonal_Prior")

  sapply(1:nrow(d), function(i) {
    val <- d[i, ]

    s <- ifelse(val["Ino"] == 1, 2, 1) # Ino = 1 represents WBFT, Ino = 2 represents EBFT

    mov <- paste0("mov_ymarrs[1, , 1, , , ", s, "]")
    start <- paste0("recdist_rs[, ", s, "]")

    r <- val["Area"]

    m <- paste0("log(", val["Index"], ")")
    cv <- val["CV"]

    paste0("calc_eqdist(", mov, ", m_start = 2, start = ", start, ")[", val["Season"], ", ", r, "] %>% log() %>% dnorm(", m, ", ", cv, ", log = TRUE)")
  })
})

# Recruitment distribution prior for stock 2 in GOM
prior_recdist <- "dnorm(p$log_recdist_rs[1, 2], 0, 1.5, log = TRUE)"

# Sample recdist prior, i.e. why sd = 1.5 seems appropriate for a mostly uninformative uniform prior away from bounds
if (FALSE) {
  nsamp <- 1e5
  nr <- 2
  set.seed(324)
  x <- rnorm(nsamp * (nr-1), 0, 1.5) %>%
    matrix(nsamp, nr-1)
  recdist <- apply(cbind(x, 0), 1, softmax)

  png("figures/recdist_prior.png", widht = 4, height = 5, units = "in", res = 400)
  par(mfrow = c(2, 1))
  hist(recdist[1, ], xlab = "Recruitment proportion WBFT in GOM", main = NULL)
  hist(recdist[2, ], xlab = "Recruitment proportion WBFT in WATL", main = NULL)
  dev.off()
}

Dmodel@prior <- c(prior_recdev, prior_dist, prior_recdist)








#### Stock configurations with biological parameters (all fixed) ----
Dstock <- new(
  "Dstock",
  m_spawn = 2,      # Spawn in season 2
  m_advanceage = 2, # Advance age classes in season 2
  SRR_s = c("BH", "BH"),
  delta_s = c(0, 0) # Spawning at start of season
)

# Length at age with Richards function and SD in length at age
# EBFT parameters reduce to von Bertalanffy equation with A2 = 999, A1 = t0 implies L1 = 0, and p = 1
growth_par <- data.frame(
  Stock = c("EBFT", "WBFT"),
  A1 = c(-0.97, 0),
  A2 = c(999, 34),
  L1 = c(0, 33),
  L2 = c(318.85, 270.6),
  K = c(0.093, 0.22),
  p = c(1, -0.12),
  a = c(3.50801e-5, 1.77054e-5),
  b = c(2.878451, 3.001252)
)
readr::write_csv(growth_par, "tables/growth_par.csv")

Dstock@len_ymas <- sapply(1:Dmodel@ns, function(s) {
  A1 <- growth_par$A1[s]
  A2 <- growth_par$A2[s]
  L1 <- growth_par$L1[s]
  L2 <- growth_par$L2[s]
  K  <- growth_par$K[s]
  p <- growth_par$p[s]

  sapply(1:Dmodel@na, function(aa) {
    sapply(1:Dmodel@nm, function(m) {
      tt <- (aa - 1) + (m - 1)/Dmodel@nm # age 0, 0.25, 0.5, ...
      ans <- L1^p + (L2^p - L1^p) * (1 - exp(-K * (tt - A1)))/(1 - exp(-K * (A2 - A1)))
      rep(ans^(1/p), Dmodel@ny)
    })
  }, simplify = "array")
}, simplify = "array")

## Override May 20, 2026:
## Use WBFT growth curve for both stocks
Dstock@len_ymas[, , , 1] <- Dstock@len_ymas[, , , 2]

Dstock@sdlen_ymas <- 0.06 * Dstock@len_ymas + 5.84

# Stock weight at age
Dstock@swt_ymas <- local({
  swt_ymas <- array(NA_real_, dim(Dstock@len_ymas))
  swt_ymas[, , , 1] <- growth_par$a[1] * Dstock@len_ymas[, , , 1] ^ growth_par$b[1]
  swt_ymas[, , , 2] <- growth_par$a[2] * Dstock@len_ymas[, , , 2] ^ growth_par$b[2]
  swt_ymas
})

# Define stock presence by area - define areas where stock can go to
Dstock@presence_rs <- matrix(FALSE, Dmodel@nr, Dmodel@ns)
Dstock@presence_rs[-1, 1] <- TRUE # EBFT can go to WATL, EATL, MED
Dstock@presence_rs[-4, 2] <- TRUE # WBFT can go to GOM, WATL, EATL

# Define areas where stocks can spawn
Dstock@natal_rs <- matrix(0, Dmodel@nr, Dmodel@ns)
Dstock@natal_rs[-1, 1] <- 1
Dstock@natal_rs[-4, 2] <- 1

#### Create two separate stock objects here:
# A: Younger maturity ogive (identical for both stocks) IN CONJUNCTION with high M
# B: Older maturity at age IN CONJUNCTION with low M and older senescence
Dstock_A <- Dstock_B <- Dstock

# Object A
Dstock_A@mat_yas <- local({
  mat_young <- rep(1, Dmodel@na)
  mat_young[1:5] <- c(0, 0, 0, 0.25, 0.5) # ages 0-4; age 5+ = 1

  array(mat_young, c(Dmodel@na, Dmodel@ny, Dmodel@ns)) %>%
    aperm(c(2, 1, 3))
})

Dstock_A@M_yas <- local({

  if (FALSE) {
    ## This is from M3, which did not include age zero
    M_high <- rep(NA_real_, Dmodel@na)
    M_high[0:14 + 1] <- c(0.38, 0.38, 0.3, 0.24, 0.2, 0.18, 0.16, 0.14, 0.13, 0.12, 0.12, 0.11, 0.11, 0.11, 0.1) # age 0-14
    M_high[15:25 + 1] <- 0.1
    M_high[seq(26, Dmodel@na)] <- 0.1

    array(M_high, c(Dmodel@na, Dmodel@ns, Dmodel@ny)) %>%
      aperm(c(3, 1, 2))
  } else {

    # From SS3 where reference M = 0.1 at age 20
    # Provided by A. Kimoto, May 20, 2026
    M_high <- array(NA_real_, c(Dmodel@ny, Dmodel@na, Dmodel@ns))

    # EBFT
    M_high[, 1:16, 1] <- matrix(
      c(0.50, 0.32, 0.27, 0.22, 0.19, 0.17, 0.15, 0.14, 0.13, 0.12, 0.12, 0.11, 0.11, 0.11, 0.11, 0.10),
      Dmodel@ny, 16, byrow = TRUE
    )

    # WBFT
    M_high[, 1:16, 2] <- matrix(
      c(0.39, 0.33, 0.29, 0.25, 0.21, 0.19, 0.17, 0.15, 0.14, 0.13, 0.12, 0.12, 0.11, 0.11, 0.11, 0.10),
      Dmodel@ny, 16, byrow = TRUE
    )

    M_high[, seq(17, Dmodel@na), ] <- 0.1

    return(M_high)
  }
})

# Object B
Dstock_B@mat_yas <- local({
  mat_E <- mat_W <- rep(1, Dmodel@na)
  mat_E[1:9] <- c(0, 0, 0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9) # EBFT age 0-8; 9+ = 1
  mat_W[1:13] <- c(0, 0, 0, 0, 0, 0, 0, 0.01, 0.04, 0.19, 0.56, 0.88, 0.98) # WBFT age 0-12; 13+ = 1

  matd_yas <- array(NA_real_, c(Dmodel@ny, Dmodel@na, Dmodel@ns))
  matd_yas[, , 1] <- matrix(mat_E, Dmodel@ny, Dmodel@na, byrow = TRUE)
  matd_yas[, , 2] <- matrix(mat_W, Dmodel@ny, Dmodel@na, byrow = TRUE)
  matd_yas
})

Dstock_B@M_yas <- local({
  M_low <- rep(NA_real_, Dmodel@na)
  M_low[0:14 + 1] <- c(0.36, 0.36, 0.27, 0.21, 0.17, 0.14, 0.12, 0.11, 0.1, 0.09, 0.09, 0.08, 0.08, 0.08, 0.08) # age 0-14
  M_low[15:25 + 1] <- 0.07
  M_low[seq(26, Dmodel@na)] <- 0.47

  array(M_low, c(Dmodel@na, Dmodel@ns, Dmodel@ny)) %>%
    aperm(c(3, 1, 2))
})











#### Fishery data - Catch (tonnes) ----
Dfishery <- new("Dfishery")
Dfishery@nf <- nrow(fleet_names)

Catch <- readxl::read_excel(xlsx_file, sheet = "Catch") %>%
  mutate(y = Year - ModelYear[1] + 1) %>%
  as.matrix()
Dfishery@Cobs_ymfr <- array(0, c(Dmodel@ny, Dmodel@nm, Dfishery@nf, Dmodel@nr))
Dfishery@Cobs_ymfr[Catch[, c("y", "Season", "Fleet", "Area")]] <- Catch[, "Catch"]

### Equilibrium catch prior to first year of model:
Dfishery@Cinit_mfr <- array(0.5 *Dfishery@Cobs_ymfr[1, , , ], c(Dmodel@nm, Dfishery@nf, Dmodel@nr))


#### Fishery data - CAL ----
len_bin <- readxl::read_excel(xlsx_file, sheet = "Length_classes")
CAL <- readxl::read_excel(xlsx_file, sheet = "Length_Comp") %>%
  mutate(y = Year - ModelYear[1] + 1) %>%
  as.matrix()
stopifnot(length(unique(CAL[, "Len_class"])) == nrow(len_bin))

Dfishery@CALobs_ymlfr <- array(0, c(Dmodel@ny, Dmodel@nm, Dmodel@nl, Dfishery@nf, Dmodel@nr))

Dfishery@CALobs_ymlfr[CAL[, c("y", "Season", "Len_class", "Fleet", "Area")]] <- CAL[, "N"]
Dfishery@CALN_ymfr <- apply(Dfishery@CALobs_ymlfr, c(1, 2, 4, 5), sum)
Dfishery@CALN_ymfr <- pmin(Dfishery@CALN_ymfr, 50)

Dfishery@fcomp_like <- "lognormal"

Dfishery@sel_f <- ifelse(fleet_names$Selectivity == "Logistic", "logistic_length", "dome_length")

stopifnot(all(apply(Dfishery@CALobs_ymlfr, 4, sum) > 0)) # Every fleet has length samples?


#### Indices of abundance: CPUE and fishery-independent indices ----
Dsurvey <- new("Dsurvey")

# From M3 trial spec document
len_cat <- data.frame(
  Name = c("US_RR_66_114", "US_RR_115_144", "US_RR_177", "US_RR_145", "US_RR_195", "US_RR_66_144"),
  Lmin = c(66, 115, 177, 145, 195, 66),
  Lmax = c(114, 144, Inf, Inf, Inf, 144)
) %>%
  mutate(bin_min = Dmodel@lmid[findInterval(Lmin, Dmodel@lmid)],
         bin_max = Dmodel@lmid[findInterval(Lmax, Dmodel@lmid)])

cpue <- readxl::read_excel(xlsx_file, sheet = "CPUE") %>%
  mutate(CV = as.numeric(CV), Index = as.numeric(Index))

cpue_use <- c("SPN_BB", "SPN_FR_BB", "MOR_SPN_TRAP", "MOR_POR_TRAP", "JPN_LL_Eatl_Med",
              "JPN_LL_NEAtl1", "JPN_LL_NEAtl2", "US_RR_66_144", "US_RR_177", "US_RR_145",
              "US_RR_195", "MEXUS_GOM_PLL", "JPN_LL_West1", "JPN_LL_West2", "JPLL_GOM",
              "CAN SWNS")
cpue %>% filter(Name %in% cpue_use, is.na(Index) | is.na(CV))

cpue_names <- summarise(
  cpue,
  n = n(),
  Fleet = unique(Fleet),
  Area = unique(Area),
  .by = Name
) %>%
  mutate(Name2 = paste0("(", c(LETTERS, letters[1:4]), ") ", Name)) %>%
  filter(Name %in% cpue_use) %>%
  mutate(i = match(Name, cpue_use)) %>%
  left_join(select(len_cat, Name, bin_min, bin_max), by = "Name") %>%
  mutate(Fleet2 = ifelse(is.na(bin_min), Fleet, paste0(Fleet, "_", bin_min, "_", bin_max)))

cpue_value <- cpue %>%
  filter(Name %in% cpue_use) %>%
  filter(!is.na(Index)) %>%
  mutate(CV = ifelse(is.na(CV), 0.3, CV)) %>% # Assume this is for MOR_POR_TRAP where CV is consistently around 0.3 for other years
  mutate(y = Year - ModelYear[1] + 1,
         i = match(Name, cpue_names$Name) |> as.numeric()) %>%
  select(!Name) %>%
  as.matrix()

index <- readxl::read_excel(xlsx_file, sheet = "Survey") %>%
  mutate(CV = as.numeric(CV), Index = as.numeric(Index)) %>%
  filter(!is.na(Index))

index_use <- c("FR_AER_SUV1", "FR_AER_SUV2", "MED_LAR_SUV", "GBYP_AER_SUV_BAR",
               "GOM_LAR_SUV", "CAN_ACO_SUV1", "CAN_ACO_SUV2")
index_names <- summarise(
  index,
  n = n(),
  Area = unique(Area),
  Stock = unique(Stock),
  Type = unique(Type),
  .by = Name
) %>%
  mutate(Name2 = paste0("(", 1:nrow(.), ") ", Name)) %>%
  filter(Name %in% index_use) %>%
  mutate(i = max(cpue_value[, "i"]) + as.numeric(match(Name, index_use)))

index_value <- index %>%
  filter(Name %in% index_use) %>%
  mutate(CV = ifelse(is.na(CV), 0.75, CV)) %>% # Assume this is for GBYP_EAR_SUV_BAR in 2024, CV in 2023 is 0.75
  mutate(y = Year - ModelYear[1] + 1,
         i = max(cpue_value[, "i"]) + as.numeric(match(Name, index_names$Name))) %>%
  select(!Name & !Type) %>%
  as.matrix()

Dsurvey@ni <- length(unique(cpue_value[, "i"])) + length(unique(index_value[, "i"]))
Dsurvey@Iobs_ymi <- Dsurvey@Isd_ymi <- array(NA_real_, c(Dmodel@ny, Dmodel@nm, Dsurvey@ni))

Dsurvey@Iobs_ymi[cpue_value[, c("y", "Season", "i")]] <- cpue_value[, "Index"]
Dsurvey@Iobs_ymi[index_value[, c("y", "Season", "i")]] <- index_value[, "Index"]

Dsurvey@Isd_ymi[cpue_value[, c("y", "Season", "i")]] <- cpue_value[, "CV"]
Dsurvey@Isd_ymi[index_value[, c("y", "Season", "i")]] <- index_value[, "CV"]

Dsurvey@unit_i <- rep("B", Dsurvey@ni) # All indices have biomass units

# Identify area and stock that each index samples. 1 = TRUE for index i in area r for stock s
cpue_samp <- array(0, c(max(cpue_names$i), nr, ns))
cpue_names_matrix <- cpue_names %>%
  select(i, Area) %>%
  as.matrix()
cpue_samp <- sapply2(1:ns, function(s) {
  cpue_samp <- array(0, c(max(cpue_names$i), nr))
  cpue_samp[cpue_names_matrix[, c("i", "Area")]] <- 1
  return(cpue_samp)
})

index_samp <- array(0, c(length(unique(index_names$i)), nr, ns))
index_names_matrix <- index_names %>%
  mutate(i = i - max(cpue_names[, "i"])) %>%
  select(!Name & !Name2 & !Type) %>%
  as.matrix()
index_samp[index_names_matrix[, c("i", "Area", "Stock")]] <- 1

Dsurvey@samp_irs <- abind::abind(cpue_samp, index_samp, along = 1)

# Selectivity of indices
# For CPUE, identify the fleet, also size range if necessary
# For stock-specific indices, identify either B and SB
Dsurvey@sel_i <- c(cpue_names$Fleet2, index_names$Type)

mutate(cpue_names, Fleet_mirror = fleet_names$FleetName[cpue_names$Fleet])

# Assume sampling at beginning of the season
Dsurvey@delta_i <- 0










#### Tag transitions ----
etag <- readr::read_csv(file.path("data", "Etag", "Etag_proportions_04.26.2026.csv")) %>%
  mutate(AgeName = ageclass_key$Name[AgeClass]) #%>%
  #summarise(N = sum(N), Nfr = sum(Nfr), .by = c(Stock, Quarter, From, To)) %>%
  #mutate(p = N/Nfr, .by = c(Stock, Quarter, From))

etag_matrix <- etag %>%
  mutate(s = ifelse(Stock == "EBFT", 1, 2),
         y = 1,
         fr = match(From, area_names$Name) %>% as.numeric(),
         to = match(To, area_names$Name) %>% as.numeric()) %>%
  select(y, s, Quarter, fr, to, N, Nfr, AgeClass, p) %>%
  as.matrix()

Dtag <- new("Dtag")

Dtag@tag_ymarrs <- array(0, c(1, Dmodel@nm, 3, Dmodel@nr, Dmodel@nr, Dmodel@ns))
Dtag@tag_ymarrs[etag_matrix[, c("y", "Quarter", "AgeClass", "fr", "to", "s")]] <- etag_matrix[, "N"]

# Data informs all years equally (constant movement with years)
Dtag@tag_yy <- matrix(1, 1, Dmodel@ny)

# Three age classes in dataset
Dtag@tag_aa <- matrix(0, 3, Dmodel@na)
Dtag@tag_aa[1, 0:4 + 1] <- Dtag@tag_aa[2, 5:8 + 1] <- Dtag@tag_aa[3, 10:Dmodel@na] <- 1

# Multinomial distribution with sample size
Dtag@tag_like <- "multinomial"
Dtag@tagN_ymars <- apply(Dtag@tag_ymarrs, c(1, 2, 3, 4, 6), sum)












#### Fishery data - Stock of origin ----
#### DECISION: Only include SOO data in the WATL and EATL!

# Prep the data object further:
# Three age classes 0-4, 5-8, 9+
Dfishery@SC_aa <- matrix(0, 3, Dmodel@na)
Dfishery@SC_aa[1, 0:4 + 1] <- Dfishery@SC_aa[2, 5:8 + 1] <- Dfishery@SC_aa[3, 10:Dmodel@na] <- 1

Dfishery@SC_like <- "logitnormal"

#### Separate dummy fleets by SOO type

# Dummy fleets: 1st row = microchemistry, 2nd row = genetics
Dfishery@SC_ff <- matrix(1, 2, Dfishery@nf)

# Otolith
SOO_otolith <- read.csv(file.path("data", "SOO", "Isotope_mixing_Proportion_Estimates_v2.csv")) %>%
  filter(Region %in% c("WATL", "EATL")) %>%
  mutate(
    EBFT = N * (1 - Prob_West),
    WBFT = N - EBFT,
    a = ifelse(grepl("0-4", fAGE), 1, ifelse(grepl("5-8", fAGE), 2, 3)),
    r = match(Region, area_names$Name),
    quarter = substr(Quarter, 2, 2) |> as.integer(),
    y = fYear - ModelYear[1] + 1
  ) %>%
  reshape2::melt(
    id.vars = c("a", "y", "quarter", "r", "SE", "N"),
    measure.vars = c("EBFT", "WBFT")
  ) %>%
  mutate(Stock = ifelse(variable == "EBFT", 1, 2), f = 1) %>%
  select(!variable) %>%
  as.matrix()
stopifnot(all(!is.na(SOO_otolith)))

# Genetic
SOO_genetic <- read.csv(file.path("data", "SOO", "Genetic_mixing_Proportion_Estimates.csv")) %>%
  filter(Region %in% c("WATL", "EATL")) %>%
  mutate(
    EBFT = N * (1 - Prob_West),
    WBFT = N - EBFT,
    f = 2,
    a = ifelse(grepl("0-4", fAGE), 1, ifelse(grepl("5-8", fAGE), 2, 3)),
    r = match(Region, area_names$Name),
    quarter = substr(Quarter, 2, 2) |> as.integer(),
    y = fYear - ModelYear[1] + 1
  ) %>%
  reshape2::melt(
    id.vars = c("a", "y", "quarter", "r", "SE", "N"),
    measure.vars = c("EBFT", "WBFT")
  ) %>%
  mutate(Stock = ifelse(variable == "EBFT", 1, 2), f = 2) %>%
  select(!variable) %>%
  as.matrix()
stopifnot(all(!is.na(SOO_genetic)))


# 3 age classes and 2 dummy fleets
Dfishery@SC_ymafrs <- array(0, c(Dmodel@ny, Dmodel@nm, 3, 2, Dmodel@nr, Dmodel@ns))
Dfishery@SC_ymafrs[SOO_otolith[, c("y", "quarter", "a", "f", "r", "Stock")]] <- SOO_otolith[, "value"]
Dfishery@SC_ymafrs[SOO_genetic[, c("y", "quarter", "a", "f", "r", "Stock")]] <- SOO_genetic[, "value"]

Dfishery@SCstdev_ymafrs <- array(NA, dim(Dfishery@SC_ymafrs))
Dfishery@SCstdev_ymafrs[SOO_otolith[, c("y", "quarter", "a", "f", "r", "Stock")]] <- SOO_otolith[, "SE"]
Dfishery@SCstdev_ymafrs[SOO_genetic[, c("y", "quarter", "a", "f", "r", "Stock")]] <- SOO_genetic[, "SE"]

















#### Labels for plotting
Dlabel <- new(
  "Dlabel",
  year = ModelYear,
  season = paste("Season", 1:4),
  age = 1:Dmodel@na - 1,
  region = area_names$Name,
  stock = c("EBFT", "WBFT"),
  fleet = fleet_names$FleetName,
  index = c(cpue_names$Name2, index_names$Name2)
)







#### Save objects
dir_save <- "model_input/05.20.2026"
if (!dir.exists(dir_save)) dir.create(dir_save)

saveRDS(Dmodel, file.path(dir_save, "Dmodel.rds"))
saveRDS(Dstock_A, file.path(dir_save, "Dstock_A.rds"))
saveRDS(Dstock_B, file.path(dir_save, "Dstock_B.rds"))
saveRDS(Dfishery, file.path(dir_save, "Dfishery.rds"))
saveRDS(Dsurvey, file.path(dir_save, "Dsurvey.rds"))
saveRDS(Dtag, file.path(dir_save, "Dtag.rds"))
saveRDS(Dlabel, file.path(dir_save, "Dlabel.rds"))
