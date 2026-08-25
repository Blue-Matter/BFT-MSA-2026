
library(r4ss)
library(tidyverse)

# SSB SS3
dirs <- c("Early", "R39")
reps <- lapply(1:length(dirs), function(i) {
  replist <- r4ss::SS_output(
    file.path("..", "BFT_2026", dirs[i]),
    verbose = FALSE
  )
})

SSB_SS3 <- lapply(1:length(dirs), function(i) {
  reps[[i]]$timeseries %>%
    filter(Era == "TIME") %>%
    select(Yr, SpawnBio) %>%
    mutate(area = ifelse(i == 1, "WATL", "EATL"),
           model = "SS3") %>%
    rename(year = Yr, S = SpawnBio)
}) %>%
  bind_rows() %>%
  select(year, area, model, S)
readr::write_csv(SSB_SS3, file = "tables/SS3_SSB.csv")

g <- ggplot(SSB_SS3, aes(year, value)) +
  geom_line() +
  facet_wrap(vars(area))
