
library(sf)
library(tidyverse)

BFT1=cbind(x=c(-80,-88,-95,-100,-100,-85,-80, -80), y=c(20,20,16.5,20,35,35,25, 20))
BFT2=cbind(x=c(-82.5,-75,-75,-65,-65,-55,-55,-70,-95,-88,-80,-80,-82.5),
          y=c(30,30,25,25,20,20,0,0,16.5,20,20,25,30))
BFT3=cbind(x=c(-70,-70,-60,-55,-55, -70), y=c(45,55,55,50,45, 45))
BFT4=cbind(x=c(-70,-55,-55,-65,-65,-75,-75,-82.5,-85,-70,-55,-55,-60,-70,-80,-100,-100,-45,-45,-30,-30,-25,-25,-70, -70),
          y=c(0,0,20,20,25,25,30,30,35,45,45,50,55,55,50,60,80,80,10,10,5,5,-50,-50, 0))
BFT5=cbind(x=c(-30,-45,-45,-30, -30), y=c(40,40,80,80, 40))
BFT6=cbind(x=c(-30,-45,-45,-30, -30), y=c(10,10,40,40, 10))
BFT7=cbind(x=c(-30,45,45,15,15,-15,-15,-30,-30), y=c(80,80,50,50,60,60,50,50,80))
BFT8=cbind(x=c(-30,-30,-15,-15,15,15,5,-5, -30), y=c(40,50,50,60,60,50,50,40, 40))
BFT9=cbind(x=c(-30,-30,-5,-5,20,20,-25,-25,-30, -30), y=c(10,40,40,30,30,-50,-50,5,5, 10))
BFT10=cbind(x=c(-5,-5,5,23,23, -5), y=c(30,40,50,50,30, 30))
BFT11=cbind(x=c(23,45,45,23, 23) ,y=c(50,50,30,30, 50))

BFT_polygon <- list(BFT1, BFT2, BFT3, BFT4, BFT5, BFT6, BFT7,
                    BFT8, BFT9, BFT10, BFT11) %>%
  lapply(list) %>%
  lapply(st_polygon)

BFT_poly <- do.call(st_sfc, BFT_polygon) %>%
  st_sf(
    stratum = c("GOM", "WATL", "WATL", "WATL", "EATL", "EATL", "EATL", "EATL", "EATL",
                "MED", "MED"),
    id = 1:11,
    crs = 4326
  )
stopifnot(all(st_is_valid(BFT_poly)))
saveRDS(BFT_poly, "data/Etag/stratum_polygon.rds")

coast <- rnaturalearth::ne_countries()
g <- ggplot(BFT_poly, aes(fill = factor(id))) +
  geom_sf() +
  geom_sf(data = coast, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  labs(fill = NULL) +
  coord_sf(xlim = c(-100, 45), ylim = c(-50, 80))



# Assign stratum for each row of geolocations
dat <- read.csv('data/Etag/BFT_geolocations_2026_03_31.csv') %>%
  arrange(group, tag, year, month, day)

dat_points <- lapply(1:nrow(dat), function(i) {
  st_point(x = c(dat$lon[i], dat$lat[i]))
})
dat_sf <- do.call(st_sfc, dat_points) %>%
  st_sf(
    geometry = .,
    dat,
    crs = 4326
  )

dat_stratum <- st_intersects(dat_sf, BFT_poly)
dat_sf$area <- sapply(dat_stratum, function(i) BFT_poly$stratum[i])

dat_sf$area <- sapply(dat_stratum, function(i) {
  if (!length(i)) return(NA)
  BFT_poly$stratum[i[length(i)]]
})

coast <- rnaturalearth::ne_countries()
g <- ggplot() +
  #geom_sf(data = BFT_poly, aes(fill = stratum)) +
  #geom_sf(data = coast) +
  geom_sf(data = dat_sf, aes(colour = area)) +
  scale_fill_viridis_d() +
  labs(fill = NULL) +
  coord_sf(xlim = c(-100, 45), ylim = c(-50, 80))

# Aggregate by transitions
tag_id <- unique(dat$tag)
tag_summary <- lapply(tag_id, function(i) {
  dat_i <- filter(dat_sf, tag == i)
  df <- data.frame()
  if (length(unique(dat_i$group)) > 1) {
    #print(i)
    dat_i <- filter(dat_i, group == dat_i$group[1])
  }

  j_start <- 1
  stratum_start <- dat_i$area[j_start] %>% as.character()
  date_start <- paste(dat_i$month[j_start], dat_i$day[j_start], dat_i$year[j_start], sep = "/") %>%
    as.POSIXct(format = "%m/%d/%Y", tz = "UTC")

  for (j in 2:nrow(dat_i)) {
    if (j >= nrow(dat_i) || dat_i$area[j] != stratum_start) {
      if (j >= nrow(dat_i)) {
        j_end <- nrow(dat_i)
      } else {
        j_end <- j-1
      }
      date_end <- paste(dat_i$month[j_end], dat_i$day[j_end], dat_i$year[j_end], sep = "/") %>%
        as.POSIXct(format = "%m/%d/%Y", tz = "UTC")

      df_new <- data.frame(
        Group = unique(dat_i$group),
        Tag_ID = unique(dat_i$tag),
        Size_at_release = unique(dat_i$size_at_tagging),
        Wt_at_release = unique(dat_i$wt_at_tagging),
        Stock_Area = stratum_start,
        Days = as.numeric(date_end - date_start) + 1,
        Start_Date = date_start,
        End_Date = date_end
      )

      df <- rbind(df, df_new)

      if (j < nrow(dat_i)) {
        j_start <- j
        stratum_start <- dat_i$area[j]
        date_start <- paste(dat_i$month[j], dat_i$day[j], dat_i$year[j], sep = "/") %>%
          as.POSIXct(format = "%m/%d/%Y", tz = "UTC")
      }
    }

  }

  return(df)
}) %>%
  bind_rows() %>%
  mutate(Entry = 1:n())

readr::write_csv(tag_summary, "data/Etag/BFT_etags_processed_forMSA_20261002.csv")
