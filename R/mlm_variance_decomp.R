#' Confidence and prediction intervals for simple slopes in a random-slope model
#'
#' In a model with a random slope for `pred`, two different questions can be
#' asked about the slope at a given moderator value \eqn{w}:
#'
#' 1. **How precisely is the average slope known?** The confidence interval
#'    for \eqn{\beta_1 + \beta_3 w} reflects only estimation uncertainty in
#'    the fixed effects (as in [mlm_probe()]).
#' 2. **What slope should be expected in a new cluster?** The prediction
#'    interval adds the residual random-slope variance \eqn{\tau_{11}}:
#'    \deqn{\hat\beta_1 + \hat\beta_3 w \pm t_{df}\sqrt{SE^2(w) + \hat\tau_{11}}.}
#'
#' Because `modx` is in the model, \eqn{\tau_{11}} is the slope variance
#' that remains *after* the moderator: a significant interaction does not mean
#' the moderator explains the between-cluster slope differences, and a large
#' \eqn{\tau_{11}} relative to the average slope means individual clusters
#' can have slopes far from it.
#'
#' @param model An `lmerMod` object with a random slope for `pred` and a
#'   two-way interaction between `pred` and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param modx.values Strategy for moderator values. See `mlm_probe()`.
#' @param at Optional numeric vector of custom moderator values.
#' @param conf.level Level for the confidence and prediction intervals.
#'   Default `0.95`.
#' @inheritParams mlm_probe
#'
#' @details
#' The prediction interval treats \eqn{\hat\tau_{11}} as known and assumes
#' normally distributed random slopes, so it is approximate and too narrow
#' when the number of clusters is small. The random slope is taken from the
#' grouping factor whose random-effects terms include `pred`.
#'
#' Versions before 0.3.0 reported `pct_random`, the ratio
#' \eqn{\tau_{11} / (\tau_{11} + SE^2)}. It mixed a population variance with
#' a sampling variance, so it grew with sample size without any change in
#' the data-generating process; it has been removed.
#'
#' @return An object of class `mlm_variance_decomp` (a list) with:
#'   * `decomp`: data frame with columns `modx_value`, `slope`, `se_fixed`
#'     (SE of the average slope), `tau11` (random-slope variance), `tau11_sd`
#'     (its square root), `se_total` (\eqn{\sqrt{SE^2 + \tau_{11}}}), `df`,
#'     `ci_lower`, `ci_upper` (confidence interval for the average slope), and
#'     `pi_lower`, `pi_upper` (prediction interval for a new cluster).
#'   * `tau11`, `tau11_sd`: the random-slope variance and SD.
#'   * `has_random_slope`: logical; does the model include a random slope for
#'     `pred`?
#'   * Metadata: `pred`, `modx`, `conf.level`, `df_method`, `grp_name`.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200), x = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 + x | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' vd <- mlm_variance_decomp(mod, pred = "x", modx = "m")
#' print(vd)
#'
#' @export
mlm_variance_decomp <- function(model,
                                pred,
                                modx,
                                modx.values = c("mean-sd", "quartiles",
                                                "tertiles", "custom"),
                                at         = NULL,
                                conf.level = 0.95,
                                df_method  = c("satterthwaite", "kenward-roger",
                                               "between", "residual"),
                                modx.level = c("auto", "cluster", "observation")) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  modx.values <- match.arg(modx.values)
  df_method   <- match.arg(df_method)
  modx.level  <- match.arg(modx.level)

  ctx       <- .infer_ctx(model, df_method)
  modx_vals <- .pick_modx_values(.modx_vector(model, modx, modx.level),
                                 modx.values = modx.values, at = at)

  # Random-slope variance from the grouping factor that carries `pred`.
  vc  <- lme4::VarCorr(model)
  grp <- NULL
  for (g in names(vc)) {
    if (pred %in% rownames(as.matrix(vc[[g]]))) { grp <- g; break }
  }
  has_random_slope <- !is.null(grp)
  if (is.null(grp)) grp <- .cluster_name(model)
  tau11 <- if (has_random_slope) as.matrix(vc[[grp]])[pred, pred] else 0

  rows <- lapply(modx_vals, function(w) {
    ss <- .simple_slope(ctx, pred, modx, w, conf.level = conf.level)
    se_total <- sqrt(ss$se^2 + tau11)
    tc <- stats::qt(1 - (1 - conf.level) / 2, df = ss$df)
    data.frame(
      modx_value = w,
      slope      = ss$slope,
      se_fixed   = ss$se,
      tau11      = tau11,
      tau11_sd   = sqrt(tau11),
      se_total   = se_total,
      df         = ss$df,
      ci_lower   = ss$ci_lower,
      ci_upper   = ss$ci_upper,
      pi_lower   = ss$slope - tc * se_total,
      pi_upper   = ss$slope + tc * se_total,
      stringsAsFactors = FALSE
    )
  })
  decomp_df <- do.call(rbind, rows)
  rownames(decomp_df) <- NULL

  structure(
    list(
      decomp           = decomp_df,
      tau11            = tau11,
      tau11_sd         = sqrt(tau11),
      has_random_slope = has_random_slope,
      pred             = pred,
      modx             = modx,
      modx.values      = modx.values,
      conf.level       = conf.level,
      df_method        = df_method,
      grp_name         = grp
    ),
    class = "mlm_variance_decomp"
  )
}

#' @export
print.mlm_variance_decomp <- function(x, digits = 3, ...) {

  cat("\n========================================\n")
  cat("  Slope Variance Decomposition\n")
  cat("========================================\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Cluster grouping:", x$grp_name, "\n")
  cat("Confidence level:", x$conf.level, "\n\n")

  if (!x$has_random_slope) {
    cat("NOTE: No random slope for '", x$pred, "' found in this model.\n",
        "Prediction intervals equal confidence intervals.\n",
        "Consider fitting (1 + ", x$pred, " | cluster) for full decomposition.\n\n",
        sep = "")
  } else {
    cat(sprintf("Random-slope SD (tau11):  %.3f\n", x$tau11_sd))
    cat(sprintf("Random-slope Var (tau11): %.3f\n\n", x$tau11))
  }

  cat("--- Per-moderator-value intervals ---\n\n")
  df <- x$decomp
  pct <- round(100 * x$conf.level)
  out <- data.frame(
    modx     = round(df$modx_value, digits),
    slope    = round(df$slope, digits),
    se       = round(df$se_fixed, digits),
    df       = round(df$df, 1),
    ci       = sprintf("[%.*f, %.*f]", digits, df$ci_lower, digits, df$ci_upper),
    pi       = sprintf("[%.*f, %.*f]", digits, df$pi_lower, digits, df$pi_upper)
  )
  names(out) <- c(x$modx, "slope", "SE", "df",
                  paste0(pct, "% CI (average slope)"),
                  paste0(pct, "% PI (new cluster)"))
  print(out, row.names = FALSE)
  if (x$has_random_slope) {
    cat("\nThe prediction interval treats tau11 as known; it is approximate\n",
        "with few clusters.\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Plot the variance decomposition of simple slopes
#'
#' Shows the simple slope at each moderator value as a point with two
#' interval layers: an inner confidence interval (fixed-effect uncertainty
#' only) and an outer prediction interval (fixed + random-slope variance).
#' The gap between the two intervals is the contribution of random-slope
#' heterogeneity.
#'
#' @param x An `mlm_variance_decomp` object.
#' @param x_label x-axis label. Defaults to moderator name.
#' @param y_label y-axis label. Defaults to "Simple slope of pred".
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.mlm_variance_decomp <- function(x,
                                     x_label = NULL,
                                     y_label = NULL,
                                     ...) {

  if (is.null(x_label)) x_label <- x$modx
  if (is.null(y_label)) y_label <- paste0("Simple slope of '", x$pred, "'")

  df <- x$decomp
  df[[x$modx]] <- df$modx_value   # for aes mapping

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = modx_value, y = slope)) +

    # Outer: prediction interval (fixed + random)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pi_lower, ymax = pi_upper),
      width = 0.06, color = "#7B2D8B", linewidth = 0.8,
      linetype = "dashed"
    ) +

    # Inner: confidence interval (fixed only)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.04, color = "#2166AC", linewidth = 1.2
    ) +

    # Slope estimate
    ggplot2::geom_point(size = 3.5, color = "#2166AC") +

    # Zero reference
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.6) +

    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = "Simple Slope Variance Decomposition",
      subtitle = paste0(
        "Inner bars: fixed-effect CI  |  ",
        "Outer bars: prediction interval (adds \u03C411 = ",
        round(x$tau11_sd, 3), ")"
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = 10)
    )

  p
}
