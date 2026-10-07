library(mlmoderator)
library(lme4)

# Shared fixtures ---------------------------------------------------------
make_mod <- function(J = 30, n = 20, seed = 1, slope_sd = 0.4) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(J), each = n)),
                  x = rnorm(J * n), w = rep(rnorm(J), each = n),
                  z = rnorm(J * n))
  d$y <- 0.3 * d$x + 0.4 * d$x * d$w + 0.2 * d$z +
    rep(rnorm(J, 0, slope_sd), each = n) * d$x +
    rep(rnorm(J), each = n) + rnorm(J * n)
  assign("fixture_data", d, envir = globalenv())
  suppressMessages(lme4::lmer(y ~ x * w + z + (1 + x | g), data = fixture_data))
}
