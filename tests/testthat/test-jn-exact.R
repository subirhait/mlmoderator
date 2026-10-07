mod <- make_mod(J = 40, seed = 3)

test_that("closed-form JN roots satisfy |t| = t_crit (constant df)", {
  jn <- mlm_jn(mod, "x", "w", df_method = "between", modx.range = c(-10, 10))
  expect_true(length(jn$jn_bounds_all) >= 1)
  for (b in jn$jn_bounds_all) {
    ss <- mlm_probe(mod, "x", "w", at = b, df_method = "between")$slopes
    expect_equal(abs(ss$t), qt(0.975, ss$df), tolerance = 1e-8)
  }
})

test_that("Satterthwaite JN boundaries satisfy the exact criterion", {
  jn <- mlm_jn(mod, "x", "w", modx.range = c(-10, 10))
  skip_if(all(is.na(jn$jn_bounds)))
  for (b in jn$jn_bounds) {
    ss <- mlm_probe(mod, "x", "w", at = b)$slopes
    expect_equal(abs(ss$t), qt(0.975, ss$df), tolerance = 1e-6)
  }
})

test_that("JN print describes the significant region", {
  jn <- mlm_jn(mod, "x", "w", df_method = "between", modx.range = c(-10, 10))
  expect_output(print(jn), "Johnson-Neyman")
  expect_s3_class(plot(jn), "ggplot")
})
