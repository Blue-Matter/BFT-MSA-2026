
# This code assigns the correct areas and processes tag tracks into single seasonal transitions
# Set up cluster for parallel computing (now 20x faster!)

# R Script for formatting PSAT data
# Updated by Q. Huynh April 2026

library(snowfall)
library(lubridate)
library(tidyverse)

# --- Matt Lauretta --- GBYP - DFO - NOAA - WWF - Unimar - FCP - CB - IEO - AZTI - UCA
dat<-read.csv("data/Etag/BFT_etags_processed_forM3_20260424.csv")


# Exceptions ---------
# Matt Lauretta erroneous fish
#cond<-dat$Tag_ID=="519703100" & dat$Stock_Area=="GOM"
#dat$Stock_Area[cond]<-"W_ATL"
#cond<-dat$Tag_ID=="2010_55308" & dat$Stock_Area=="GOM"
#dat$Stock_Area[cond]<-"CAR"


#### Format dates
dat <- dat[dat$Days > 0, ]
dat$Start_Date <- as.POSIXct(dat$Start_Date, format = "%m/%d/%Y", tz = "UTC")
dat$End_Date <- as.POSIXct(dat$End_Date, format = "%m/%d/%Y", tz = "UTC")

cond <- is.na(dat$Start_Date)

if (any(cond)) {
  datstr<-dat
  dat$Start_Date[cond]<-as.POSIXct(datstr$Start_Date[cond], format = "%m/%d/%Y", tz = "UTC")
  dat$End_Date[cond]<-as.POSIXct(datstr$End_Date[cond], format = "%m/%d/%Y", tz = "UTC")
  rm(datstr)
}

#### These are inputs now hardwired in this code (they were formally slots of the OMI object)
np = 2  # nstocks
nr = 4  # nareas
na = 35 # nages
nma = 3 # nmovementageclasses
ns = 4 # nseasons (impSOO.R only)

len_age = c(53.37763,  76.95330,  98.43530, 118.00956, 135.84549, 152.09748, 166.90619, 180.39979, 192.69507, 203.89845, 214.10689,
223.40876, 231.88456, 239.60766, 246.64489, 253.05718, 258.90001, 264.22396, 269.07510, 273.49544, 277.52322, 281.19330,
284.53746, 287.58464, 290.36121, 292.89120, 295.19651, 297.29709, 299.21113, 300.95519, 302.54436, 303.99241, 305.31186,
306.51414, 307.60964)

#wt_age = c(3.289878, 9.429172, 19.153491, 32.282927, 48.409102, 67.018108, 87.566675, 109.528128, 132.418566, 155.810104, 179.335775, 202.689213, 225.621200, 247.934508, 269.477948,
#290.140243, 309.844091, 328.540637, 346.204473, 362.829203, 378.423589, 393.008229, 406.612752, 419.273460, 431.031371, 441.930621, 452.017158, 461.337707, 469.938943, 477.866857,
#485.166266, 491.880454, 498.050918, 503.717195, 508.916752)

mov_ac = c(rep(1,4),rep(2,4),rep(3,na-8))

areanams = c("GOM", "WATL", "EATL", "MED")

#### Cohort slice: size at release to age
# Should be okay as we will aggregate into larger age classes
getlb <- function(x, vec) {
  vec <- c(0, vec)
  max((1:length(vec))[x>vec])
}

age <- sapply(dat$Size_at_release, getlb, vec = len_age)

if (any(is.na(age))) stop("There are NAs when cohort slicing Size_at_release to age")
#agew <- sapply(dat$Weight_kg, getlb, vec = wt_age)
#age[is.na(age)]<-agew[is.na(age)]

#### Aggregate integer ages to age classes
dat$ageclass <- mov_ac[age]

#### Check age class and start/end dates
dat0 <- subset(dat,!(is.na(dat$ageclass)|is.na(dat$Start_Date)|is.na(dat$End_Date)))

#### Assign stock of origin to individual tags: WBFT if they spent more days in GOM than MED
#Areas <- c(match("MED", areanams), match("GOM", areanams))
#Natal<- All[, c("TagID", "Area")]
#Natal <- Natal[Natal$Area %in% c("MED", "GOM"), ]
#NatalIDs <- data.frame(ID = Natal$TagID, Stock = ifelse(Natal$Area == "MED", "EBFT", "WBFT"))
#NatalIDs <- aggregate(NatalIDs$Stock, by = list(NatalIDs$ID), max)
#names(NatalIDs)<-c("TagID","Stock")

NatalIDs <- dat0 %>%
  filter(Stock_Area %in% c("MED", "GOM")) %>%
  summarise(Days = sum(Days), .by = c(Tag_ID, Stock_Area)) %>%
  pivot_wider(names_from = Stock_Area, values_from = Days, values_fill = 0) %>%
  mutate(Stock = ifelse(GOM > MED, "WBFT", "EBFT"))

#filter(NatalIDs, GOM == MED)
#NotNatal <- filter(dat0, !Tag_ID %in% NatalIDs$Tag_ID)


#### Expands a record (within area) into a daily set of records
tagexpand <- function(r, dat0) {
  nd <- as.integer((dat0$End_Date[r] - dat0$Start_Date[r]) + 1)
  Date <- dat0$Start_Date[r] + lubridate::days(seq(0, nd-1))
  Year <- as.numeric(format(Date, "%Y"))
  Quarter <- ceiling(as.numeric(format(Date, "%m"))/3)

  data.frame(
    TagID=as.character(dat0$Tag_ID[r]),
    Area=dat0$Stock_Area[r],
    Date=Date,
    Year=Year,
    Quarter=Quarter,
    AgeC=dat0$ageclass[r],
    stringsAsFactors = FALSE
  )
}

# Init parallel processing
sfInit(cpus = parallel::detectCores() - 2, parallel = TRUE)

temp_L<-sfLapply(1:nrow(dat0),tagexpand,dat0=dat0)
All<-bind_rows(temp_L,.id = "column_label")

# As per M3 format population, subyear, time duration (quarters) til capture, from area, to area, N
#stk<-array(rep(1:np,each=nr),c(nr,np))
#ar<-array(rep(1:nr,np),c(nr,np))
#ExSpawn<-cbind(c(1,2),c(4,1)) # Exclusive spawning areas (for Stock ID of tags) each row is a stock, second column is the exclusive natal area

TagIDs<-unique(All$TagID)
nTags<-length(TagIDs)

# Assume all AZTI tags are Eastern origin
#defosEAST<-unique(dat$Tag_ID[dat$Stock_Area==6 & dat$Investigator=="AZTI"]) # Dr Haritz Arrizabalaga stipulated 'certain Eastern fish'

defosEAST <- dat0 %>% filter(Group == "AZTI") %>% pull(Tag_ID) %>% unique()
#filter(dat0, Tag_ID %in% defosEAST)

#### Simplify tracks to quarterly transitions

# Find the area in which a tag spent the most days in a quarterly time step
dat1 <- summarise(All, ndays = n(), .by = c(TagID, Year, Quarter, AgeC, Area)) %>%
  summarise(Area = Area[which.max(ndays)[1]], .by = c(TagID, Year, Quarter, AgeC)) %>%
  arrange(Year, Quarter) %>%
  rename(AgeClass = AgeC)

# Convert dat1 to transition summary
simp_tracks <- function(i, dat1, TagIDs, byyear = TRUE, defosEAST, NatalIDs) {

  Trk <- subset(dat1, dat1$TagID == TagIDs[i])
  nT <- nrow(Trk)

  if (nT>1) {
    pop <- NatalIDs$Stock[match(TagIDs[i], NatalIDs$Tag_ID)]
    if (is.na(pop) && TagIDs[i] %in% defosEAST) pop <- "EBFT" # if EATL and AZTI (Haritz' stipulation)

    transitions <- data.frame(From = Trk$Area[seq(1, nT-1)], To = Trk$Area[seq(2, nT)])

    outtrack <- data.frame(Stock = rep(pop, nT - 1)) %>%
      cbind(Trk[-nT, c("AgeClass", "Year", "Quarter")]) %>%
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

Track_L <- sfLapply(1:nTags, simp_tracks, dat1 = dat1, TagIDs = TagIDs, byyear = TRUE,
                    defosEAST = defosEAST, NatalIDs = NatalIDs)
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
complete <- apply(Tracks, 1, function(i) all(!is.na(i)))
print(paste(round(sum(complete)/nrow(Tracks)*100,1),"% of tags are of known stock of origin. Removing ",nrow(Tracks)-sum(complete)," tracks."))
Tracks <- Tracks[complete, ] # remove any line with unknown stock, subyear, or area

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

readr::write_csv(PSAT, "data/Etag/Etag_proportions_04.26.2026.csv")


