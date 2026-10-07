mod <- make_mod(J = 12, n = 15, seed = 5)

test_that("LOCO refits reproduce manual refits with the same REML setting", {
  s <- mlm_sensitivity(mod, "x", "w", df_method = "between")
  expect_s3_class(s, "mlm_sensitivity")
  expect_equal(nrow(s$loco), 12L)
  expect_equal(s$summary$n_failed, 0L)
  fd <- lme4::getData(mod)
  d1 <- fd[fd$g != "1", ]
  f1 <- suppressMessages(lme4::lmer(y ~ x * w + z + (1 + x | g), data = d1))
  expect_true(lme4::isREML(f1))
  expect_equal(s$loco$b_int[s$loco$cluster == "1"],
               unname(lme4::fixef(f1)["x:w"]), tolerance = 1e-6)
  expect_equal(s$loco$dfbeta,
               (s$loco$b_int - s$observed$b_int) / s$observed$se_int)
  expect_equal(s$cutoff, 2 / sqrt(12))
  expect_output(print(s), "Leave-one-cluster-out")
  expect_s3_class(plot(s), "ggplot")
})

test_that("LOCO works with transformed variables in the formula", {
  fd <- lme4::getData(mod)
  fd$ypos <- exp(fd$y / 5)
  m2 <- suppressMessages(lme4::lmer(log(ypos) ~ x * w + (1 | g), data = fd))
  s <- mlm_sensitivity(m2, "x", "w", df_method = "between")
  expect_equal(s$summary$n_failed, 0L)
})

test_that("defunct arguments warn or error", {
  expect_warning(mlm_sensitivity(mod, "x", "w", df_method = "between",
                                 icc_range = c(0, .3)), "removed")
  expect_error(mlm_sensitivity(mod, "x", "w", loco = FALSE), "no longer")
})
