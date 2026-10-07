mod <- make_mod()

test_that("mlm_surface returns a ggplot", {
  expect_s3_class(mlm_surface(mod, "x", "w", grid = 20), "ggplot")
  expect_s3_class(mlm_surface(mod, "x", "w", grid = 20, fill = FALSE,
                              probe_lines = FALSE), "ggplot")
})

test_that("mlm_summary bundles interaction, slopes, and JN", {
  s <- mlm_summary(mod, "x", "w", df_method = "between")
  expect_s3_class(s, "mlm_summary")
  expect_equal(s$interaction$term, "x:w")
  expect_equal(nrow(s$simple_slopes), 3L)
  expect_s3_class(s$jn, "mlm_jn")
  expect_output(print(s), "Interaction Term")
  expect_null(mlm_summary(mod, "x", "w", jn = FALSE,
                          df_method = "between")$jn)
})

test_that("plot confidence bands account for all fixed effects", {
  p <- mlm_plot(mod, "x", "w", df_method = "between")
  g <- p$data
  X <- mlmoderator:::.fixed_X(mod, g)
  se <- sqrt(rowSums((X %*% as.matrix(vcov(mod))) * X))
  tc <- qt(0.975, 28)
  expect_equal(unname(g$ci_hi - g$fitted), unname(tc * se), tolerance = 1e-8)
})
