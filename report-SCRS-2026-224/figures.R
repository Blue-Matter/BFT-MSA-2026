
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
