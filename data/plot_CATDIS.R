
library(tidyverse)

cdis <- readr::read_csv("data/cdis5024-BFT.csv")

sumry <- cdis %>%
  summarise(Catch_t = sum(Catch_t), .by = c(Stock, Decade, yLat5ctoid, xLon5ctoid))

borders <- rnaturalearth::ne_countries()
g <- ggplot(sumry, aes(xLon5ctoid, yLat5ctoid)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_point(aes(colour = log(Catch_t))) +
  facet_wrap(vars(Decade)) +
  scale_colour_viridis_c() +
  coord_sf(xlim = c(-100, 50), ylim = c(-50, 70))

g <- sumry %>%
  filter(Stock == "ATW") %>%
  ggplot(aes(xLon5ctoid, yLat5ctoid)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_point(aes(colour = log(Catch_t))) +
  facet_wrap(vars(Decade)) +
  scale_colour_viridis_c() +
  coord_sf(xlim = c(-100, -20), ylim = c(-50, 70))

catch_ATW <- cdis %>%
  filter(Stock == "ATW") %>%
  mutate(Fleet = ifelse(PartyName %in% c("UNITED STATES", "JAPAN", "CANADA"), PartyName,
                        ifelse(GearGrp == "LL", "LLOTH",
                               ifelse(GearGrp == "PS", "PS", "OTHER")))) %>%
  summarise(Catch_t = sum(Catch_t), .by = c(Fleet, Decade, Stock, yLat5ctoid, xLon5ctoid))

g <- catch_ATW %>%
  ggplot(aes(xLon5ctoid, yLat5ctoid)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_tile(width = 5, height = 5, aes(fill = log(Catch_t))) +
  facet_grid(vars(Decade), vars(Fleet)) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-100, -20), ylim = c(-50, 70))

catch_ATW_all <- cdis %>%
  filter(Stock == "ATW") %>%
  mutate(Fleet = ifelse(GearGrp == "LL" & PartyName != "JAPAN", "LL_OTHER",
                        ifelse(GearGrp == "PS", "PS",
                               ifelse(PartyName %in% c("UNITED STATES", "JAPAN", "CANADA"), PartyName, "OTHER")))) %>%
  summarise(Catch_t = sum(Catch_t), .by = c(Fleet, yLat5ctoid, xLon5ctoid))

g <- catch_ATW_all %>%
  ggplot(aes(xLon5ctoid, yLat5ctoid)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_tile(width = 5, height = 5, aes(fill = log(Catch_t))) +
  facet_wrap(vars(Fleet)) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-100, -20), ylim = c(-50, 70))
ggsave("figures/data/CATDIS_ATW.png", g, height = 6, width = 6)


d <- data.frame(
  Decade = seq(1950, 2020, 10),
  Score = c("1950-1969", "1970-1989", "1990-2009", "2010-") |> rep(each = 2)
)
catch_ATW_all <- cdis %>%
  left_join(d, by = "Decade") %>%
  filter(Stock == "ATW") %>%
  mutate(Fleet = ifelse(GearGrp == "LL" & PartyName != "JAPAN", "LL_OTHER",
                        ifelse(GearGrp == "PS", "PS",
                               ifelse(PartyName %in% c("UNITED STATES", "JAPAN", "CANADA"), PartyName, "OTHER")))) %>%
  summarise(Catch_t = sum(Catch_t), .by = c(Fleet, Score, yLat5ctoid, xLon5ctoid))

g <- catch_ATW_all %>%
  ggplot(aes(xLon5ctoid, yLat5ctoid)) +
  geom_sf(data = borders, inherit.aes = FALSE) +
  geom_tile(width = 5, height = 5, aes(fill = log(Catch_t))) +
  facet_grid(vars(Score), vars(Fleet)) +
  scale_fill_viridis_c() +
  coord_sf(xlim = c(-100, -20), ylim = c(-50, 70)) +
  labs(x = "Longitude", y = "Latitude")
ggsave("figures/data/CATDIS_ATW_score.png", g, height = 8, width = 9)
