mod <- make_mod()

test_that("prediction interval adds tau11 to the fixed-effect variance", {
  vd <- mlm_variance_decomp(mod, "x", "w", df_method = "between")
  tau11 <- as.matrix(lme4::VarCorr(mod)$g)["x", "x"]
  expect_true(vd$has_random_slope)
  expect_equal(vd$tau11, tau11)
  d <- vd$decomp
  expect_equal(d$se_total, sqrt(d$se_fixed^2 + tau11))
  tc <- qt(0.975, d$df)
  expect_equal(d$pi_upper - d$slope, tc * d$se_total)
  expect_true(all(d$pi_upper - d$pi_lower > d$ci_upper - d$ci_lower))
  expect_false("pct_random" %in% names(d))
  expect_output(print(vd), "PI \\(new cluster\\)")
  expect_s3_class(plot(vd), "ggplot")
})

test_that("without a random slope, PI equals CI", {
  m0 <- suppressMessages(lme4::lmer(y ~ x * w + (1 | g), data = fixture_data))
  vd <- mlm_variance_decomp(m0, "x", "w", df_method = "between")
  expect_false(vd$has_random_slope)
  expect_equal(vd$decomp$pi_lower, vd$decomp$ci_lower)
})
