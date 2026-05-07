

library(tidyverse)

dir_save <- "model_input/04.30.2026"

Dstock_A <- readRDS(file.path(dir_save, "Dstock_A.rds"))
Dstock_B <- readRDS(file.path(dir_save, "Dstock_B.rds"))

# Length at age
len_age <- Dstock_A@len_ymas[1, , , ] %>%
  reshape2::melt() %>%
  rename(Season = Var1, Int_Age = Var2, Length = value, Stock = Var3) %>%
  mutate(Age = (Int_Age - 1) + 0.25 * (Season - 1), Stock = ifelse(Stock == 1, "EBFT", "WBFT"))

sdlen_age <- Dstock_A@sdlen_ymas[1, , , ] %>%
  reshape2::melt() %>%
  rename(Season = Var1, Int_Age = Var2, SD = value, Stock = Var3) %>%
  mutate(Age = (Int_Age - 1) + 0.25 * (Season - 1), Stock = ifelse(Stock == 1, "EBFT", "WBFT"))

g <- left_join(len_age, sdlen_age) %>%
  ggplot(aes(Age, Length, fill = Stock, colour = Stock)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Length - 1.96 * SD,
                  ymax = Length + 1.96 * SD),
              linetype = NA,
              alpha = 0.25) +
  expand_limits(y = 0) +
  labs(y = "Length")
ggsave("figures/bio/length-age.png", g, width = 5, height = 3)

# Weight at age
wt_age <- Dstock_A@swt_ymas[1, , , ] %>%
  reshape2::melt() %>%
  rename(Season = Var1, Int_Age = Var2, Weight = value, Stock = Var3) %>%
  mutate(Age = (Int_Age - 1) + 0.25 * (Season - 1), Stock = ifelse(Stock == 1, "EBFT", "WBFT"))

g1 <- left_join(len_age, wt_age) %>%
  ggplot(aes(Length, Weight, fill = Stock, colour = Stock)) +
  geom_line() +
  expand_limits(y = 0)

g2 <- wt_age %>%
  ggplot(aes(Age, Weight, fill = Stock, colour = Stock)) +
  geom_line() +
  expand_limits(y = 0) +
  labs(y = "Weight")

g <- ggpubr::ggarrange(g1, g2, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave("figures/bio/weight-age.png", g, width = 6, height = 3)

# Maturity at age
Mat_age_A <- Dstock_A@mat_yas[1, , ] %>%
  reshape2::melt() %>%
  rename(Age = Var1, Maturity = value, Stock = Var2) %>%
  mutate(Age = Age - 1, Stock = ifelse(Stock == 1, "EBFT", "WBFT"), Type = "(A) Early maturity, high M")

Mat_age_B <- Dstock_B@mat_yas[1, , ] %>%
  reshape2::melt() %>%
  rename(Age = Var1, Maturity = value, Stock = Var2) %>%
  mutate(Age = Age - 1, Stock = ifelse(Stock == 1, "EBFT", "WBFT"), Type = "(B) Late maturity, low M")


M_A <- Dstock_A@M_yas[1, , ] %>%
  reshape2::melt() %>%
  rename(Age = Var1, M = value, Stock = Var2) %>%
  mutate(Age = Age - 1, Stock = ifelse(Stock == 1, "EBFT", "WBFT"), Type = "(A) Early maturity, high M")

M_B <- Dstock_B@M_yas[1, , ] %>%
  reshape2::melt() %>%
  rename(Age = Var1, M = value, Stock = Var2) %>%
  mutate(Age = Age - 1, Stock = ifelse(Stock == 1, "EBFT", "WBFT"), Type = "(B) Late maturity, low M")

g <- left_join(
  rbind(Mat_age_A, Mat_age_B),
  rbind(M_A, M_B)
) %>%
  reshape2::melt(id.vars = c("Age", "Stock", "Type")) %>%
  ggplot(aes(Age, value, colour = Stock)) +
  facet_grid(vars(variable), vars(Type), scales = "free_y", switch = "y") +
  geom_line() +
  geom_point() +
  labs(y = NULL) +
  expand_limits(y = 0) +
  theme(strip.placement = "outside", strip.background = element_blank())
ggsave("figures/bio/maturity-M.png", g, width = 6, height = 4)



