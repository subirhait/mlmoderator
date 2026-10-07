mod <- make_mod()

test_that("simple slopes match hand computation (between/residual df)", {
  fe <- lme4::fixef(mod); V <- as.matrix(vcov(mod))
  w <- 0.7
  est <- fe[["x"]] + w * fe[["x:w"]]
  se  <- sqrt(V["x", "x"] + w^2 * V["x:w", "x:w"] + 2 * w * V["x", "x:w"])
  pr <- mlm_probe(mod, "x", "w", at = w, df_method = "between")$slopes
  expect_equal(pr$slope, est, tolerance = 1e-10)
  expect_equal(pr$se, se, tolerance = 1e-10)
  # J = 30 clusters, one cluster-level predictor (w) -> 30 - 1 - 1
  expect_equal(pr$df, 28)
  pr2 <- mlm_probe(mod, "x", "w", at = w, df_method = "residual")$slopes
  expect_equal(pr2$df, nrow(mod@frame) - length(fe))
})

test_that("Satterthwaite results match lmerTest and emmeans", {
  mt <- lmerTest::as_lmerModLmerTest(mod)
  ref <- summary(mt)$coefficients["x:w", ]
  s <- mlm_summary(mod, "x", "w", jn = FALSE)$interaction
  expect_equal(s$df, unname(ref["df"]), tolerance = 1e-6)
  expect_equal(s$p, unname(ref["Pr(>|t|)"]), tolerance = 1e-6)
  skip_if_not_installed("emmeans")
  at <- c(-1, 0, 1.2)
  em <- summary(emmeans::emtrends(mod, ~ w, var = "x", at = list(w = at),
                                  lmer.df = "satterthwaite"))
  pr <- mlm_probe(mod, "x", "w", at = at)$slopes
  expect_equal(pr$slope, em$x.trend, tolerance = 1e-6)
  expect_equal(pr$se, em$SE, tolerance = 1e-6)
  expect_equal(pr$df, em$df, tolerance = 1e-3)
})

test_that("Kenward-Roger matches lmerTest when pbkrtest is available", {
  skip_if_not_installed("pbkrtest")
  mt <- lmerTest::as_lmerModLmerTest(mod)
  ref <- summary(mt, ddf = "Kenward-Roger")$coefficients["x:w", ]
  s <- mlm_summary(mod, "x", "w", jn = FALSE, df_method = "kenward-roger")$interaction
  expect_equal(s$df, unname(ref["df"]), tolerance = 1e-4)
  expect_equal(s$se, unname(ref["Std. Error"]), tolerance = 1e-6)
})

test_that("cluster-level moderator values are computed over clusters", {
  fd <- lme4::getData(mod)
  w_cl <- tapply(fd$w, fd$g, `[`, 1)
  pr <- mlm_probe(mod, "x", "w", df_method = "between")
  expect_equal(pr$slopes$modx_value,
               mean(w_cl) + c(-1, 0, 1) * sd(w_cl), tolerance = 1e-10)
  pr_obs <- mlm_probe(mod, "x", "w", df_method = "between",
                      modx.level = "observation")
  expect_equal(pr_obs$slopes$modx_value[2], mean(fd$w))
  expect_error(mlm_probe(mod, "w", "x", modx.level = "cluster"), "varies within")
})

test_that("df_method is validated", {
  expect_error(mlm_probe(mod, "x", "w", df_method = "wrong"))
})
