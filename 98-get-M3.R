
library(ABTMSE) # version 8.2.4
library(parallel)

ABTMSE::loadABT()

ABTMSE::Design$all_lnams

OM_subset <- ABTMSE::Design$LNames_Ref %>%
  mutate(n = 1:n()) %>%
  filter(grepl("all years", Var1),
         grepl("A", Var2),
         grepl("H", Var4)) %>%
  mutate(label = as.character(Var3) %>% strsplit(": ") %>% sapply(getElement, 1))

OM <- lapply(OM_subset$n, function(i) get(paste0("OM_", i, "d")))

cl <- parallel::makeCluster(length(OM))
parallel::clusterEvalQ(cl, library(ABTMSE))
parallel::clusterEvalQ(cl, ABTMSE::loadABT())
MSElist <- parallel::parLapply(
  cl,
  OM,
  function(i) {
    i@interval <- 200L
    new("MSE", OM = i)
  }
)
parallel::stopCluster(cl)
saveRDS(MSElist, file = "M3/ABTMSE_list.rds")

SSB <- lapply(1:length(MSElist), function(i) {
  M <- MSElist[[i]]
  M@SSB[1, 1, , 1:M@nyears] %>%
    structure(
      dimnames = list(stock = ifelse(M@Snames == "East", "EBFT", "WBFT"),
                      year = seq(1, M@nyears) + 2018 - M@nyears)
    ) %>%
    reshape2::melt(value.name = "S") %>%
    mutate(model = paste("M3:", OM_subset$label[i]))
}) %>%
  bind_rows()
readr::write_csv(SSB, file = "tables/M3_SSB.csv")

