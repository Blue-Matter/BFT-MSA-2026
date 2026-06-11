

library(multiSA)


model_name <- c(
  "reference1_06.03.2026.rds",
  "reference2_06.03.2026.rds",
  "Wprior1_06.03.2026.rds",
  "Wprior2_06.03.2026.rds"
)

#for (j in 1:length(model_name)) {
for (j in 3:4) {
  i <- model_name[j]
  fit <- readRDS(file.path("model_output", paste0("fit_", i)))

  tictoc::tic()
  fits <- do_jitter(fit, n = 10, amount = 0.25, seed = 23, cores = 3, use_fitted = FALSE, return_models = TRUE, do_sd = FALSE)
  tictoc::toc()

  saveRDS(fits, file = file.path("jitter", paste0("jitter_", i)))

  # Need time to compute max gradient
  like <- do.call(rbind, lapply(fits, get_likelihood_components))
  saveRDS(like, file = file.path("jitter", paste0("jitter_loglik_", i)))
}

if (FALSE) {
  i <- model_name[4]
  like <- readRDS(file = file.path("jitter", paste0("jitter_loglik_", i)))

  g <- like %>%
    mutate(Jitter = 1:n()) %>%
    reshape2::melt(id.vars = "Jitter") %>%
    filter(!variable %in% c("loglike_Cinit_mfr", "conv", "fn")) %>%
    mutate(value = ifelse(grepl("log", variable), -1 * value, value)) %>%
    #mutate(value = value - min(value), .by = variable) %>%
    ggplot(aes(Jitter, value)) +
    facet_wrap(vars(variable), scales = "free_y") +
    geom_point() +
    geom_line() +
    labs(x = "Jitter Run")

}
