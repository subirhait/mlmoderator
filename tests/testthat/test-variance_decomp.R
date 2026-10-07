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
  m0 <- suppressMessages(lme4::lmer(y ~ x * w + (1 | g), data = lme4::getData(mod)))
  vd <- mlm_variance_decomp(m0, "x", "w", df_method = "between")
  expect_false(vd$has_random_slope)
  expect_equal(vd$decomp$pi_lower, vd$decomp$ci_lower)
})

test_that("plot subtitle reports the random-slope SD and variance correctly", {
  vd <- mlm_variance_decomp(mod, "x", "w", df_method = "between")
  p <- plot(vd)
  # ggplot2 >= 3.5.2 resolves labels lazily; get_labs() is the stable accessor
  sub <- if (exists("get_labs", asNamespace("ggplot2"))) {
    ggplot2::get_labs(p)$subtitle
  } else {
    p$labels$subtitle
  }
  expect_match(sub, paste0("SD = ", format(round(vd$tau11_sd, 3), nsmall = 3)),
               fixed = TRUE)
  expect_match(sub, paste0("variance = ", format(round(vd$tau11, 3), nsmall = 3)),
               fixed = TRUE)
  expect_false(grepl("[^ -~\n]", sub))   # ASCII only
})
