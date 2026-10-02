
# This code assigns the correct areas and processes tag tracks into single seasonal transitions
# Set up cluster for parallel computing (now 20x faster!)

# R Script for formatting PSAT data
# Updated by Q. Huynh October 2026

library(snowfall)
library(lubridate)
library(tidyverse)

# --- Matt Lauretta --- GBYP - DFO - NOAA - WWF - Unimar - FCP - CB - IEO - AZTI - UCA
dat <- read.csv("data/Etag/BFT_etags_processed_forMSA_20261002.csv")

stopifnot(all(dat$Days > 0))
stopifnot(all(!is.na(dat$Start_Date)))

#### These are inputs now hardwired in this code (they were formally slots of the OMI object)
np = 2  # nstocks
nr = 4  # nareas
na = 35 # nages
nma = 3 # nmovementageclasses
ns = 4 # nseasons (impSOO.R only)

#len_age = c(53.37763,  76.95330,  98.43530, 118.00956, 135.84549, 152.09748, 166.90619, 180.39979, 192.69507, 203.89845, 214.10689,
#223.40876, 231.88456, 239.60766, 246.64489, 253.05718, 258.90001, 264.22396, 269.07510, 273.49544, 277.52322, 281.19330,
#284.53746, 287.58464, 290.36121, 292.89120, 295.19651, 297.29709, 299.21113, 300.95519, 302.54436, 303.99241, 305.31186,
#306.51414, 307.60964)

#wt_age = c(3.289878, 9.429172, 19.153491, 32.282927, 48.409102, 67.018108, 87.566675, 109.528128, 132.418566, 155.810104, 179.335775, 202.689213, 225.621200, 247.934508, 269.477948,
#290.140243, 309.844091, 328.540637, 346.204473, 362.829203, 378.423589, 393.008229, 406.612752, 419.273460, 431.031371, 441.930621, 452.017158, 461.337707, 469.938943, 477.866857,
#485.166266, 491.880454, 498.050918, 503.717195, 508.916752)

Growth <- data.frame(
  age = seq(1, na)
)
Growth$length <- local({
  A1 <- 0
  A2 <- 34
  L1 <- 33
  L2 <- 270.6
  K <- 0.22
  p <- -0.12

  ans <- L1^p + (L2^p - L1^p) * (1 - exp(-K * (Growth$age - A1)))/(1 - exp(-K * (A2 - A1)))
  ans^(1/p)
})
Growth$weightW <- 1.77e-5 * Growth$length^3.001
Growth$weightE <- 3.51e-5 * Growth$length^2.878
Growth$age_class <- ifelse(Growth$age < 8, 1, 2)

# At what age and weight are fish 175 cm?
approx(Growth$length, y = Growth$age, xout = 175)
approx(Growth$length, y = Growth$weightW, xout = 175)
approx(Growth$length, y = Growth$weightE, xout = 175)


areanams = c("GOM", "WATL", "EATL", "MED")

# Fit length-weight relationship in tags
dat$Wt_at_release[dat$Group == "IFREMER" & dat$Tag_ID == 200069] <- 67
dat_wt <- dat %>%
  filter(!is.na(Size_at_release), !is.na(Wt_at_release)) %>%
  summarise(wt = unique(as.numeric(Wt_at_release)), size = unique(Size_at_release), .by = c(Group, Tag_ID))
plot(wt ~ size, data = dat_wt)

fit <- nls(wt ~ a * size^b, data = dat_wt, start = list(b = 3, a = 1e-3))
points(dat_wt$size, predict(fit), col = 2)

# Fish > 175 cm at about 94 kg
filter(dat_wt, size > 175) %>%
  arrange(wt) %>%
  pull(wt) %>% min()
abline(h = 94)
abline(v = 175)

#### Impute length from weight (for plotting only)
dat$Size_impute <- ifelse(!is.na(dat$Size_at_release), dat$Size_at_release,
                          (as.numeric(dat$Wt_at_release)/coef(fit)["a"])^(1/coef(fit)["b"]))


#### Aggregate to age classes, note there are tags without size at release (remove them later)
dat$ageclass <- ifelse(dat$Size_impute >= 175, 2, 1)

#### Assign stock of origin to individual tags: WBFT if they spent more days in GOM than MED
# Also assume all AZTI tags are Eastern-origin
NatalIDs <- dat %>%
  filter(Stock_Area %in% c("MED", "GOM") | Group == "AZTI") %>%
  summarise(Days = sum(Days), .by = c(Group, Tag_ID, Stock_Area)) %>%
  pivot_wider(names_from = Stock_Area, values_from = Days, values_fill = 0) %>%
  mutate(Stock = ifelse(GOM > 0, "WBFT", ifelse(MED > 0 | Group == "AZTI", "EBFT", NA)))

filter(NatalIDs, is.na(Stock))
filter(NatalIDs, Group == "AZTI", MED == 0)
filter(NatalIDs, GOM == MED)
filter(NatalIDs, MED > 0 & GOM > 0)
#NotNatal <- filter(dat0, !Tag_ID %in% NatalIDs$Tag_ID)

#### Assign stock of origin to entire database
dat_all <- left_join(dat, select(NatalIDs, Group, Tag_ID, Stock), by = c("Group", "Tag_ID")) %>%
  mutate(Release_Area = Stock_Area[which.min(Entry)], .by = Tag_ID)
filter(dat_all, Stock_Area == "MED", is.na(Stock))
filter(dat_all, Stock_Area == "GOM", is.na(Stock))

rel_area <- dat_all %>%
  summarise(Release_Area = unique(Release_Area),
            Size_impute = min(Size_impute), .by = c(Tag_ID, Group))

tracks_day <- readr::read_csv("data/Etag/BFT_geolocations_2026_03_31.csv") %>%
  left_join(rel_area, by = c("group" = "Group", "tag" = "Tag_ID")) %>%
  left_join(select(NatalIDs, Group, Tag_ID, Stock), by = c("group" = "Group", "tag" = "Tag_ID")) %>%
  filter(!is.na(Release_Area)) %>% # These are duplicates
  mutate(date = paste(month, day, year, sep = "/") |> as.Date(format = "%m/%d/%Y")) %>%
  mutate(days_at_liberty = date - min(date), .by = tag)

borders <- rnaturalearth::ne_countries()
BFT_poly <- readRDS("data/Etag/stratum_polygon.rds")

# Tracks by stock of origin
g <- tracks_day %>%
  ggplot(aes(lon, lat, group = tag)) +
  geom_sf(data = BFT_poly, inherit.aes = FALSE, alpha = 0.25, aes(fill = stratum)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_path(alpha = 0.25, linewidth = 0.25) +
  coord_sf(xlim = c(-100, 50), ylim = c(10, 70)) +
  facet_wrap(vars(Stock), ncol = 2) +
  guides(colour = "none") +
  labs(x = "Longitude", y = "Latitude", fill = NULL) +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom")
ggsave("figures/data/Etag_tracks.png", g, height = 4, width = 6)

# Tracks by stock of origin and release area
g <- tracks_day %>%
  ggplot(aes(lon, lat, group = tag)) +
  geom_sf(data = BFT_poly, inherit.aes = FALSE, alpha = 0.25, aes(fill = stratum)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_path(alpha = 1, linewidth = 0.25) +
  coord_sf(xlim = c(-100, 50), ylim = c(10, 70)) +
  facet_grid(vars(paste("Release:", Release_Area)), vars(Stock)) +
  guides(colour = "none") +
  labs(x = "Longitude", y = "Latitude", fill = NULL) +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom")
ggsave("figures/data/Etag_tracks2.png", g, height = 8, width = 6)

#### Identify tags by stock of origin, release area (+/- 175 cm) & size at release
g <- tracks_day %>%
  mutate(Size_Bin = cut(Size_impute, breaks = c(0, 175, Inf), labels = c("<175 cm", "175+ cm"))) %>%
  ggplot(aes(lon, lat, group = tag)) +
  geom_path(linewidth = 0.25, aes(colour = Release_Area)) +
  geom_sf(data = BFT_poly, inherit.aes = FALSE, alpha = 0.25) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  coord_sf(xlim = c(-100, 40), ylim = c(10, 70), expand = FALSE) +
  facet_grid(vars(Size_Bin), vars(Stock)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  labs(x = "Longitude", y = "Latitude", fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/data/Etag_tracks_size.png", g, height = 6, width = 8)

# Histogram of time at liberty
t_liberty <- tracks_day %>%
  mutate(Size_Bin = cut(Size_impute, breaks = c(0, 175, Inf), labels = c("<175 cm", "175+ cm"))) %>%
  summarise(days_at_liberty = max(days_at_liberty), .by = c(tag, Size_Bin, Stock))

t_label <- summarise(t_liberty, med = median(days_at_liberty), .by = c(Stock, Size_Bin))

g <- t_liberty %>%
  ggplot(aes(days_at_liberty)) +
  geom_histogram() +
  geom_label(data = t_label, aes(label = med), x = Inf, y = Inf, hjust = 'inward', vjust = 'inward') +
  facet_grid(vars(Size_Bin), vars(Stock), scales = "free_y") +
  labs(x = "Days at liberty")
ggsave("figures/data/Etag_DAL_size.png", g, height = 6, width = 6)


#### Identify tags by stock of origin, release area (25 cm bins) & size at release
g <- tracks_day %>%
  mutate(Size_Bin = cut(Size_impute, breaks = seq(50, 350, 25),
                        labels = c("50-75 cm", "75-100 cm", "100-125 cm", "125-150 cm", "150-175 cm", "175-200 cm", "200-225 cm", "225-250 cm",
                                   "250-275 cm", "275-300 cm", "300+ cm", "300+ cm"))) %>%
  ggplot(aes(lon, lat, group = tag)) +
  geom_path(linewidth = 0.25, aes(colour = Release_Area)) +
  geom_sf(data = BFT_poly, inherit.aes = FALSE, alpha = 0.25) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  coord_sf(xlim = c(-100, 40), ylim = c(10, 70), expand = FALSE) +
  facet_grid(vars(Size_Bin), vars(Stock)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  labs(x = "Longitude", y = "Latitude", fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/data/Etag_tracks_size_25cm.png", g, height = 15, width = 10)

t_liberty <- tracks_day %>%
  mutate(Size_Bin = cut(Size_impute, breaks = seq(50, 350, 25),
                        labels = c("50-75 cm", "75-100 cm", "100-125 cm", "125-150 cm", "150-175 cm", "175-200 cm", "200-225 cm", "225-250 cm",
                                   "250-275 cm", "275-300 cm", "300+ cm", "300+ cm"))) %>%
  summarise(days_at_liberty = max(days_at_liberty), .by = c(tag, Size_Bin, Stock))
g <- t_liberty %>%
  ggplot(aes(days_at_liberty)) +
  geom_histogram() +
  facet_grid(vars(Size_Bin), vars(Stock)) +
  labs(x = "Days at liberty")
ggsave("figures/data/Etag_DAL_size_25cm.png", g, height = 15, width = 10)


# Identify whether NA's could be identified
filter(dat_all, Stock == "EBFT", all(Stock_Area %in% "WATL"), .by = Tag_ID) # Every Eastern-origin fish has visited areas other than WATL
potential_WBFT <- filter(dat_all, is.na(Stock), all(Stock_Area %in% "WATL"), .by = Tag_ID)
g <- tracks_day %>%
  mutate(Size_Bin = cut(size_at_tagging, breaks = c(0, 175, Inf), labels = c("<175 cm", "175+ cm"))) %>%
  mutate(Stock2 = ifelse(tag %in% potential_WBFT$Tag_ID, "WBFT", "EBFT")) %>%
  filter(is.na(Stock) | Stock == "EBFT") %>%
  ggplot(aes(lon, lat, group = tag)) +
  geom_path(linewidth = 0.25, aes(colour = Stock2)) +
  geom_sf(data = BFT_poly, inherit.aes = FALSE, alpha = 0.25) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  coord_sf(xlim = c(-100, 40), ylim = c(10, 70), expand = FALSE) +
  facet_grid(vars(Release_Area), vars(Stock)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  labs(x = "Longitude", y = "Latitude", fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/data/Etag_tracks_large.png", g, height = 6, width = 8)

#### Expands a record (within area) into a daily set of records
tagexpand <- function(r, dat0) {
  nd <- dat0$Days[r]
  Date <- as.Date(dat0$Start_Date[r]) + lubridate::days(seq(0, nd-1))
  Year <- as.numeric(format(Date, "%Y"))
  Quarter <- ceiling(as.numeric(format(Date, "%m"))/3)

  data.frame(
    Group = dat0$Group[r],
    Tag_ID = as.character(dat0$Tag_ID[r]),
    Area = dat0$Stock_Area[r],
    Date = Date,
    Year = Year,
    Quarter = Quarter,
    AgeClass = dat0$ageclass[r],
    Stock = dat0$Stock[r],
    stringsAsFactors = FALSE
  )
}

# Init parallel processing
sfInit(cpus = parallel::detectCores() - 2, parallel = TRUE)

dat0 <- filter(dat_all, !is.na(Size_impute)) # Exclude tags without size
temp_L <- sfLapply(1:nrow(dat0), tagexpand, dat0 = dat0)
All <- bind_rows(temp_L, .id = "column_label")

sfStop()

# As per M3 format population, subyear, time duration (quarters) til capture, from area, to area, N
#stk<-array(rep(1:np,each=nr),c(nr,np))
#ar<-array(rep(1:nr,np),c(nr,np))
#ExSpawn<-cbind(c(1,2),c(4,1)) # Exclusive spawning areas (for Stock ID of tags) each row is a stock, second column is the exclusive natal area

#### Simplify tracks to quarterly transitions

# Find the area in which a tag spent the most days in a quarterly time step, exclude tags without stock of origin
dat1 <- summarise(All, ndays = n(), .by = c(Group, Tag_ID, Year, Quarter, AgeClass, Area, Stock)) %>%
  summarise(Area = Area[which.max(ndays)[1]], .by = c(Group, Tag_ID, Year, Quarter, AgeClass, Stock)) %>%
  arrange(Year, Quarter)

# Convert dat1 to transition summary
simp_tracks <- function(i, dat1, TagIDs, byyear = TRUE) {
  Trk <- filter(dat1, dat1$Tag_ID == TagIDs[i])
  nT <- nrow(Trk)

  if (nT>1) {
    transitions <- data.frame(From = Trk$Area[seq(1, nT-1)], To = Trk$Area[seq(2, nT)])

    outtrack <- Trk[-nT, c("AgeClass", "Year", "Quarter", "Stock")] %>%
      cbind(transitions) %>%
      mutate(TagNo = i)

    if (!byyear) {
      outtrack <- outtrack %>% select(!Year)
    }

  } else {
    outtrack <- data.frame()
  }

  return(outtrack)
}

TagIDs <- unique(dat1$Tag_ID)
nTags <- length(TagIDs)

Track_L <- sfLapply(1:nTags, simp_tracks, dat1 = dat1, TagIDs = TagIDs, byyear = TRUE)
Tracks_byyear <- bind_rows(Track_L)
Tracks <- select(Tracks_byyear, !Year) # Remove the year

Impute <- FALSE
if (Impute) {
  # This does the movement pattern matching of Carruthers (2017) SCRS/2016/205:
  # https://www.iccat.int/Documents/CVSP/CV073_2017/n_7/CV073072552.pdf
  # to assign stock of origin to tags that did not travel to a natal area
  # NOT used in bluefin MSE

  # Prior probability of SOO based on tag transition
  # Transitions from GOM are WBFT
  # Transitions from MED are EBFT
  # Otherwise, equal prior probability
  pw<-1/nr
  PriorSOO<-array(0,c(np,ns,nr,nr))
  PriorSOO[1,,,]<-rep(c(0,rep(pw,nr-1)),each=ns*nr)
  PriorSOO[2,,,]<-rep(c(rep(pw,nr-1),0),each=ns*nr)

  # Tracks split into those with SOO and those without
  Tracks_ImputeSOO <- filter(Tracks, is.na(Stock))

  # Tracks with SOO
  Tracks_Skip <- filter(Tracks, !is.na(Stock))
  TSOO <- summarise(Tracks_Skip, N = n(), .by = c(Stock, Quarter, From, To))

  # Calculate movement matrix based on tracks with SOO
  movSOO<-Priormov<-LikeSOO<-array(0,c(np,ns,nr,nr))
  movSOO[] <- reshape2::acast(TSOO, list("Stock", "Quarter", "From", "To"), value.var = "N")
  movSOO <- movSOO/array(apply(movSOO,1:3,sum, na.rm = TRUE),dim(movSOO))
  movSOO[is.na(movSOO)]<-0

  Priormov<-(movSOO+PriorSOO)/array(apply(movSOO+PriorSOO,1:3,sum, na.rm = TRUE),dim(movSOO))
  Priormov[is.na(Priormov)]<-0

  # Get equilibrium distribution implied in Priormov
  conv<-function(relsize=c(1,1),Priormov,ny=100,ns=4){
    recmov<-array(NA,dim(Priormov))
    nr<-dim(Priormov)[4]
    for(p in 1:np){
      vec<-rep(relsize[p]/nr,nr)
      for(y in 1:ny){
        for(s in 1:ns){
          recmov[p,s,,]<-vec*Priormov[p,s,,]
          vec<-apply(recmov[p,s,,],2,sum)
        }
      }
    }
    recmov
  }
  recmov<-conv(c(1,1),Priormov,ny=20)

  LHD1<-(recmov[1,,,])/(recmov[1,,,]+recmov[2,,,])
  LHD1[is.na(LHD1)] <- 0
  LHD1[, , areanams == "MED"] <- 1  # MED prob eastern is 1
  LHD1[, , areanams == "GOM"] <- 0  # GOM prob eastern is 0

  LHD2<-recmov[2,,,]/(recmov[1,,,]+recmov[2,,,])
  LHD2[is.na(LHD2)] <- 0
  LHD2[, , areanams == "MED"]<-0  # #MED prob western is 0
  LHD2[, , areanams == "GOM"]<-1   # GOM prob western is 1

  Imptagnos <- Tracks_ImputeSOO %>% pull(TagNo) %>% unique()
  ratio <- rep(NA_real_, length(Imptagnos))

  for (tt in 1:length(Imptagnos)) {

    temp <- filter(Tracks_ImputeSOO, TagNo == Imptagnos[tt])
    ind <- temp %>%
      select(Quarter, From, To) %>%
      mutate(From = match(From, areanams), To = match(To, areanams)) %>%
      as.matrix()

    prob1 <- prod(LHD1[ind])
    prob2 <- prod(LHD2[ind])

    ratio[tt] <- prob1/prob2
    if (ratio[tt] == Inf || is.na(ratio[tt])) ratio[tt] <- 1

    stk <- NA

    if (ratio[tt] > 2) stk <- "EBFT"      # ratio of > 2 is assigned eastern
    if (ratio[tt] < 0.5) stk <- "WBFT"    # ratio of < 0.5 is assigned western
    Tracks_ImputeSOO[Tracks_ImputeSOO$TagNo == Imptagnos[tt], "Stock"] <- stk
    print(paste(Imptagnos[tt], ratio[tt], stk , sep=" - "))

  }

  Tracks <- rbind(
    Tracks_ImputeSOO,
    Tracks_Skip
  )

  # Compare transitions from imputed SOO tags and known SOO tags
  g <- rbind(
    Tracks_ImputeSOO |> mutate(Impute = TRUE),
    Tracks_Skip |> mutate(Impute = FALSE)
  ) %>%
    summarise(N = n(), .by = c(Stock, AgeClass, Quarter, From, To, Impute)) %>%
    arrange(Stock, AgeClass, Quarter, From, To) %>%
    mutate(Nfr = sum(N), .by = c(Stock, AgeClass, Quarter, From, Impute)) %>%
    mutate(p = N/Nfr) %>%
    mutate(To_i = match(To, areanams), From_i = paste("From:", From)) %>%
    filter(Stock == "EBFT") %>%
    ggplot(aes(To_i, p, shape = Impute, linetype = Impute,
               colour = factor(AgeClass))) +
    facet_grid(vars(Quarter), vars(From_i)) +
    geom_line() +
    geom_point() +
    labs(x = "To", y = "Proportion") +
    scale_linetype_manual(values = 1:2) +
    scale_shape_manual(values = c(16, 1)) +
    scale_x_continuous(labels = areanams, breaks = 1:nr)

}

# remove tags of uncertain stock of origin (do not enter a natal area)
tags_SOO <- unique(Tracks$TagNo[!is.na(Tracks$Stock)]) |> length()
tracks_SOO <- sum(!is.na(Tracks$Stock))
complete <- apply(Tracks, 1, function(i) all(!is.na(i)))

print(paste0(tracks_SOO, " (", 100 * round(tracks_SOO/nrow(Tracks), 1), "%) tag transitions are of known stock of origin"))
print(paste0(tags_SOO, " (", 100 * round(tags_SOO/length(unique(dat_all$Tag_ID)), 1), "%) tag are of known stock of origin"))
Tracks <- Tracks[!is.na(Tracks$Stock), ] # remove any line with unknown stock, subyear, or area

# Aggregate tracks into total numbers of tags (N) for any unique transition
# Calculate Nfr: number of tags exiting area and proportion of tags that moved to each region
PSAT <- summarise(Tracks, N = n(), .by = c(Stock, AgeClass, Quarter, From, To)) %>%
  arrange(Stock, AgeClass, Quarter, From, To) %>%
  mutate(Nfr = sum(N), .by = c(Stock, AgeClass, Quarter, From)) %>%
  mutate(p = N/Nfr)

#### NOTE: M3 calculates movement at beginning of time step. Therefore, the tag transitions should be
# assigned to the next quarter. Meanwhile, MSA calculates movement at the end of the time step (leave as-is).
M3 <- FALSE
if (M3) {
  PSAT <- mutate(PSAT, Quarter = ifelse(Quarter < 4), Quarter + 1, 1)
}

readr::write_csv(PSAT, "data/Etag/Etag_proportions_10.02.2026.csv")



