
library(tidyverse)

# Logit prior ----
png("figures/logit-prior.png", height = 3, width = 6, units = "in", res = 400)
set.seed(234)
x <- rnorm(1e5, 0, 1.5)
y <- plogis(x)
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
hist(x, main = NULL, ylab = "Pr(x)")
hist(y, main = NULL, ylab = "Pr(y)", xlab = expression(y == frac(1, 1 + exp(x))))
dev.off()

# Stock mixing proportions ----
SOO <- rbind(
  readr::read_csv("data/SOO/Empirical_Profile_Stock_Predictions_JPN.csv") %>%
    select(Year, Fleet, Predicted_Value, CV, Lower_95, Upper_95) %>%
    rename(P_West_Mean = Predicted_Value) %>%
    mutate(Type = "Empirical Profile"),
  readr::read_csv("data/SOO/P_West_Year_Fleet_Marginalized.csv") %>%
    select(Year, Fleet, P_West_Mean, CV, Lower_95, Upper_95) %>%
    mutate(Type = "Marginalized")
)

N <- readr::read_csv("data/SOO/P_West_Year_Fleet_Marginalized.csv") %>%
  select(Year, Fleet, N_Obs)

left_join(SOO, N)

g <- left_join(SOO, N) %>%
  ggplot(aes(Year, P_West_Mean)) +
  geom_line(linewidth = 0.5, aes(colour = Fleet)) +
  geom_linerange(linewidth = 0.25, aes(colour = Fleet, ymin = Lower_95, ymax = Upper_95)) +
  geom_point(shape = 21, alpha = 0.75, aes(fill = Fleet, size = N_Obs)) +
  facet_wrap(vars(Type), ncol = 2) +
  theme(legend.position = "bottom") +
  labs(y = "Predicted P(West)", size = "Sample Size (N)")
ggsave("figures/data/SOO_WATL.png", g, height = 4, width = 8)
