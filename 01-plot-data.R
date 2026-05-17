
library(readxl)
library(tidyverse)
library(maps)
library(mapdata)

## Master file to organize most model settings
xlsx_file <- file.path("data", "1_M3_data", "ICCAT_MSA_Data_2026Apr_v1.xlsx")

## Spatial definitions ----
## 7 areas from M3
# 1 - GOM
# 2 - WATL
# 3 - GSL
# 4 - SATL
# 5 - NATL
# 6 - EATL
# 7 - MED

## 4 areas for MSA
# 1 - GOM
# 2 - WATL + GSL
# 3 - SATL + NATL + EATL
# 4 - MED

# Taken from BFT_Stock_Area_Assign.r

#DEFINED STOCK AREA X (LON) AND Y (LAT) BOUNDARIES
BFT1=list(x=c(-80,-88,-95,-100,-100,-85,-80), y=c(20,20,16.5,20,35,35,25))
BFT2=list(x=c(-82.5,-75,-75,-65,-65,-55,-55,-70,-95,-88,-80,-80,-82.5),
          y=c(30,30,25,25,20,20,0,0,16.5,20,20,25,30))
BFT3=list(x=c(-70,-70,-60,-55,-55), y=c(45,55,55,50,45))
BFT4=list(x=c(-70,-55,-55,-65,-65,-75,-75,-82.5,-85,-70,-55,-55,-60,-70,-80,-100,-100,-45,-45,-30,-30,-25,-25,-70),
          y=c(0,0,20,20,25,25,30,30,35,45,45,50,55,55,50,60,80,80,10,10,5,5,-50,-50))
BFT5=list(x=c(-30,-45,-45,-30), y=c(40,40,80,80))
BFT6=list(x=c(-30,-45,-45,-30), y=c(10,10,40,40))
BFT7=list(x=c(-30,45,45,15,15,-15,-15,-30,-30), y=c(80,80,50,50,60,60,50,50,80))
BFT8=list(x=c(-30,-30,-15,-15,15,15,5,-5), y=c(40,50,50,60,60,50,50,40))
BFT9=list(x=c(-30,-30,-5,-5,20,20,-25,-25,-30), y=c(10,40,40,30,30,-50,-50,5,5))
BFT10=list(x=c(-5,-5,5,23,23), y=c(30,40,50,50,30))
BFT11=list(x=c(23,45,45,23) ,y=c(50,50,30,30))

#STOCK AREAS PLOTTED ON MAP
png("figures/data/areas.png", height = 6, width = 6, res = 400, units = "in")

maps::map('worldHires',col=c('gray'),fill=T,xlim=c(-100,45),ylim=c(-50,80))
axis(1,at=seq(-100,45,5))
axis(2,at=seq(-50,80,5))
mtext('Longitude',1,line=3)
mtext('Latitude',2,line=3)
polygon(BFT1,border=1,lwd=2)
text("GOM",x=-90,y=25,font=2,col=2)
polygon(BFT2,border=1,lwd=2)
text("WATL",x=-70,y=15,font=2,col=2)
polygon(BFT3,border=1,lwd=2)
text("WATL",x=-63,y=49,font=2,col=2)
polygon(BFT4,border=1,lwd=2)
text("WATL",x=-60,y=30,font=2,col=2)
polygon(BFT5,border=1,lwd=2)
text("EATL",x=-37.5,y=60,font=2,col=2)
polygon(BFT6,border=1,lwd=2)
text("EATL",x=-37.5,y=25,font=2,col=2)
polygon(BFT7,border=1,lwd=2)
text("EATL",x=-5,y=70,font=2,col=2)
polygon(BFT8,border=1,lwd=2)
text("EATL",x=-10,y=45,font=2,col=2)
polygon(BFT9,border=1,lwd=2)
text("EATL",x=-5,y=0,font=2,col=2)
polygon(BFT10,border=1,lwd=2)
text("MED",x=10,y=40,font=2,col=2)
polygon(BFT11,border=1,lwd=2)
text("MED",x=35,y=40,font=2,col=2)
dev.off()

# Original 7 area map
png("figures/data/areas_7strata.png", height = 6, width = 6, res = 400, units = "in")

maps::map('worldHires',col=c('gray'),fill=T,xlim=c(-100,45),ylim=c(-50,80))
axis(1,at=seq(-100,45,5))
axis(2,at=seq(-50,80,5))
mtext('Longitude',1,line=3)
mtext('Latitude',2,line=3)
polygon(BFT1,border=1,lwd=2)
text("GOM",x=-90,y=25,font=2,col=2)
polygon(BFT2,border=1,lwd=2)
text("WATL",x=-70,y=15,font=2,col=2)
polygon(BFT3,border=1,lwd=2)
text("GSL",x=-63,y=49,font=2,col=2)
polygon(BFT4,border=1,lwd=2)
text("W. ATL",x=-60,y=30,font=2,col=2)
polygon(BFT5,border=1,lwd=2)
text("N. ATL",x=-37.5,y=60,font=2,col=2)
polygon(BFT6,border=1,lwd=2)
text("S. ATL",x=-37.5,y=25,font=2,col=2)
polygon(BFT7,border=1,lwd=2)
text("N. ATL",x=-5,y=70,font=2,col=2)
polygon(BFT8,border=1,lwd=2)
text("E. ATL",x=-10,y=45,font=2,col=2)
polygon(BFT9,border=1,lwd=2)
text("S. ATL",x=-5,y=0,font=2,col=2)
polygon(BFT10,border=1,lwd=2)
text("MED",x=10,y=40,font=2,col=2)
polygon(BFT11,border=1,lwd=2)
text("MED",x=35,y=40,font=2,col=2)
dev.off()


### Area names
area_names <- readxl::read_excel(xlsx_file, sheet = "Areas") %>%
  select(Area, Name) %>%
  filter(1:nrow(.) %in% seq(1, nrow(.), 2))

### Fleet names
fleet_names <- readxl::read_excel(xlsx_file, sheet = "Fleets") %>%
  mutate(FleetName = paste0(Number, " - ", Code))
















### Catch
Catch <- readxl::read_excel(xlsx_file, sheet = "Catch") %>%
  mutate(Area = factor(area_names$Name[Area], area_names$Name),
         Fleet = factor(fleet_names$FleetName[Fleet], fleet_names$FleetName))
range(Catch$Area)
range(Catch$Season)
range(Catch$Year)

# Annual catch
Cannual <- Catch %>%
  summarise(Catch = sum(Catch), .by = c(Year, Area, Fleet))

g <- ggplot(Cannual, aes(Year, Catch, fill = Fleet)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.25) +
  facet_wrap(vars(Area)) +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  labs(y = "Annual Catch") +
  guides(fill = guide_legend(ncol = 4, title = "Fleet"))
ggsave("figures/data/Cobs_annual.png", g, height = 6, width = 6)

g@facet$params$free$y <- TRUE
ggsave("figures/data/Cobs_annual2.png", g, height = 6, width = 6)

# Seasonal catch
g <- ggplot(Catch, aes(Year, Catch, fill = Fleet)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.1) +
  #geom_point() +
  facet_grid(vars(Area), vars(paste("Season", Season)),
             scales = "free_y") +
  labs(x = "Year") +
  theme(panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "Seasonal Catch") +
  guides(fill = guide_legend(ncol = 4, title = NULL))
ggsave("figures/data/Cobs_seasonal.png", g, height = 6, width = 6)

# Spool-up catch at age (1950-1964, but time series already has catch during this time period)
Hist_Catch <- readxl::read_excel(xlsx_file, sheet = "Hist_Catch")

#g <- Hist_Catch %>%
#  mutate(Year = Year + 0.25 * (Season - 1)) %>%
#  mutate(Area = factor(area_names$Name[Area], area_names$Name)) %>%
#  arrange(Year) %>%
#  filter(Age %in% seq(5, 40, 5)) %>%
#  mutate(Age = paste("Age", Age) |> factor(paste("Age", seq(5, 40, 5)))) %>%
#  ggplot(aes(Year, Catch)) +
#  geom_line() +
#  geom_point(aes(colour = factor(Season))) +
#  labs(colour = "Season") +
#  facet_grid(vars(Age), vars(Area), scales = "free_y")

g <- Hist_Catch %>%
  mutate(Year = Year + 0.25 * (Season - 1)) %>%
  mutate(Area = factor(area_names$Name[Area], area_names$Name)) %>%
  arrange(Year) %>%
  #filter(Age %in% seq(5, 40, 5)) %>%
  #mutate(Age = paste("Age", Age) |> factor(paste("Age", seq(5, 40, 5)))) %>%
  ggplot(aes(Year, Catch, fill = Age)) +
  geom_col(width = 0.25, colour = "grey40", linewidth = 0.1) +
  facet_wrap(vars(Area), scales = "free_y") +
  labs(y = "Seasonal Catch") +
  scale_fill_distiller()
ggsave("figures/data/Cobs_spoolup.png", g, height = 4, width = 6)


## CAL ----
len_bin <- readxl::read_excel(xlsx_file, sheet = "Length_classes")
CAL <- readxl::read_excel(xlsx_file, sheet = "Length_Comp") %>%
  mutate(Bin = len_bin$LengthClass[Len_class]) %>%
  mutate(Area = factor(area_names$Name[Area], area_names$Name),
         Fleet = factor(fleet_names$FleetName[Fleet], fleet_names$FleetName)) %>%
  mutate(p = N/sum(N), .by = c(Year, Season, Area, Fleet))

# Mean length
g <- CAL %>%
  summarise(mlen = weighted.mean(Bin, N), .by = c(Year, Season, Fleet, Area)) %>%
  mutate(Year = Year + 0.25 * (Season - 1)) %>%
  arrange(Year) %>%
  ggplot(aes(Year, mlen, shape = factor(Season), group = Area, colour = Area)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 1) +
  facet_wrap(vars(Fleet), ncol = 3) +
  scale_shape_manual(values = c(1, 2, 4, 16)) +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  labs(x = "Year", y = "Mean length (cm)", shape = "Quarter") +
  coord_cartesian(ylim = c(0, 300)) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
ggsave("figures/data/ML.png", g, height = 8, width = 6)

# Aggregate length comp
g <- CAL %>%
  summarise(N = sum(N), .by = c(Bin, Area, Fleet)) %>%
  mutate(p = N/sum(N), .by = c(Area, Fleet)) %>%
  ggplot(aes(Bin, p, colour = Area)) +
  geom_line() +
  geom_point() +
  labs(x = "Length bin (cm)", y = "Proportion") +
  facet_wrap(vars(Fleet),
             ncol = 3,
             scales = "free_y") +
  theme(legend.position = "bottom")
ggsave("figures/data/CAL_aggregate.png", g, height = 8, width = 6)

# Fishery CPUE
cpue <- readxl::read_excel(xlsx_file, sheet = "CPUE") %>%
  mutate(CV = as.numeric(CV), Index = as.numeric(Index))

cpue_names <- summarise(cpue, n = n(), .by = Name) %>%
  mutate(Name2 = paste0("(", c(LETTERS, letters[1:4]), ") ", Name))

g <- cpue %>%
  mutate(Area = factor(area_names$Name[Area], area_names$Name)) %>%
  mutate(Year = Year + 0.25 * (Season - 1)) %>%
  left_join(cpue_names, by = "Name") %>%
  mutate(Name2 = factor(Name2, cpue_names$Name2)) %>%
  ggplot(aes(Year, Index, colour = Area, group = Area)) +
  geom_line(linewidth = 0.1, linetype = 3, color = "grey40") +
  geom_linerange(linewidth = 0.25, aes(ymin = exp(log(Index) - 2 * CV), ymax = exp(log(Index) + 2*CV))) +
  geom_point(size = 0.25) +
  facet_wrap(vars(Name2), ncol = 4, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Fishery CPUE", colour = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/data/CPUE.png", g, height = 9, width = 8)

# Fishery-independent index
index <- readxl::read_excel(xlsx_file, sheet = "Survey") %>%
  mutate(CV = as.numeric(CV), Index = as.numeric(Index))

index_names <- summarise(index, n = n(), .by = Name) %>%
  mutate(Name2 = paste0("(", 1:nrow(.), ") ", Name))

g <- index %>%
  mutate(Area = factor(area_names$Name[Area], area_names$Name)) %>%
  mutate(Year = Year + 0.25 * (Season - 1),
         Stock = ifelse(Stock == 1, "EBFT", "WBFT")) %>%
  mutate(lwr = exp(log(Index) - 2 * CV), upr = exp(log(Index) + 2*CV)) %>%
  left_join(index_names, by = "Name") %>%
  ggplot(aes(Year, Index, colour = Area, shape = Stock)) +
  geom_line(linewidth = 0.1, linetype = 3, color = "grey40") +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  geom_point(size = 1) +
  facet_wrap(vars(Name2), scales = "free_y") +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 21)) +
  labs(x = "Year", y = "Stock-specific index", colour = NULL, shape = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/data/FI_Index.png", g, height = 4, width = 6)






# Stock of origin
SOO <- rbind(
  read.csv(file.path("data", "SOO", "Isotope_mixing_Proportion_Estimates_v2.csv")) |> mutate(Source = "Otolith"),
  read.csv(file.path("data", "SOO", "Genetic_mixing_Proportion_Estimates.csv")) |> mutate(Source = "Genetic")
) |>
  mutate(Season = substr(Quarter, 2, 2) |> as.numeric()) %>%
  mutate(Year = fYear + 0.25 * (Season - 1)) %>%
  arrange(Year)

g <- SOO %>%
  mutate(Area = factor(Region, area_names$Name)) %>%
  ggplot(aes(Year, Prob_West, colour = Source, group = Source, shape = Quarter)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 1) +
  facet_grid(vars(paste("Age:", fAGE)), vars(Area)) +
  scale_shape_manual(values = c(1, 2, 4, 16)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(1970, 2026), expand = FALSE) +
  labs(x = "Year", y = "Probability WBFT") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/data/SOO.png", g, height = 5, width = 8)

# Etag
ageclass_key <- readxl::read_excel(xlsx_file, sheet = "Age_classes") %>%
  mutate(Name = ifelse(Class == 1, "0-4", ifelse(Class == 2, "5-8", "9+"))) %>%
  summarise(n = n(), .by = c(Class, Name))
etag <- readr::read_csv(file.path("data", "Etag", "Etag_proportions_04.26.2026.csv")) %>%
  mutate(AgeClass = ageclass_key$Name[AgeClass])

g <- etag %>%
  mutate(To = match(To, area_names$Name)) %>%
  mutate(From_i = paste("From:", From) |> factor(paste("From:", area_names$Name))) %>%
  ggplot(aes(To, p, colour = Stock, shape = factor(AgeClass), linetype = factor(AgeClass))) +
  facet_grid(vars(paste("Season", Quarter)), vars(From_i)) +
  geom_line() +
  geom_point(aes(size = Nfr)) +
  labs(x = "Destination", y = "Proportion", size = "N", shape = "Age", linetype = "Age") +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_linetype_manual(values = 3:1) +
  scale_shape_manual(values = c(1, 4, 16)) +
  scale_x_continuous(labels = area_names$Name, breaks = 1:nrow(area_names)) +
  guides(colour = guide_legend(ncol = 1),
         linetype = guide_legend(ncol = 1),
         size = guide_legend(ncol = 1))
ggsave("figures/data/Etag.png", g, height = 6, width = 6)
