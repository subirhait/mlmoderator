#' Probe simple slopes from a multilevel interaction
#'
#' Computes simple slopes of a focal predictor (`pred`) at selected values of a
#' moderator (`modx`) from a two-level mixed-effects model fitted with
#' [lme4::lmer()]. Returns estimates, standard errors, *t*-values, *p*-values,
#' and confidence intervals in a tidy data frame.
#'
#' @param model An `lmerMod` object containing a two-way interaction between
#'   `pred` and `modx` in the fixed-effects structure.
#' @param pred Character scalar. Name of the focal predictor variable.
#' @param modx Character scalar. Name of the moderator variable.
#' @param modx.values Strategy for selecting moderator values. One of:
#'   * `"mean-sd"` (default): mean --- 1 SD, mean, mean + 1 SD.
#'   * `"quartiles"`: 25th, 50th, 75th percentiles.
#'   * `"tertiles"`: 33rd and 67th percentiles.
#'   * `"custom"`: use values supplied via `at`.
#' @param at Numeric vector of custom moderator values. Used when
#'   `modx.values = "custom"`, or to override any strategy.
#' @param conf.level Confidence level for intervals. Default `0.95`.
#' @param df_method Denominator degrees of freedom for tests and intervals:
#'   * `"satterthwaite"` (default): Satterthwaite approximation via
#'     \pkg{lmerTest}.
#'   * `"kenward-roger"`: Kenward-Roger approximation (requires
#'     \pkg{pbkrtest}); also uses the Kenward-Roger adjusted covariance, so
#'     standard errors can differ slightly.
#'   * `"between"`: \eqn{J - q - 1}, where \eqn{J} is the number of clusters
#'     and \eqn{q} the number of cluster-level fixed effects (Snijders and
#'     Bosker, 2012). A simple, conservative choice for cross-level terms.
#'   * `"residual"`: \eqn{N - p}, the rule used by mlmoderator before 0.3.0.
#'     It treats level-1 observations as independent and is anti-conservative
#'     for cross-level interactions, especially with few clusters; it is kept
#'     only for reproducing earlier results.
#' @param modx.level How moderator values are summarised when choosing probe
#'   points. `"auto"` (default) uses one value per cluster when `modx` is
#'   constant within clusters (a cluster-level moderator) and all
#'   observations otherwise; `"cluster"` and `"observation"` force either.
#'   Ignored when `at` is supplied.
#'
#' @references
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin*, 2(6), 110--114.
#'
#' Snijders, T. A. B., & Bosker, R. J. (2012). *Multilevel analysis* (2nd ed.).
#' Sage.
#'
#' Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017). lmerTest
#' package: Tests in linear mixed effects models. *Journal of Statistical
#' Software*, 82(13), 1--26. \doi{10.18637/jss.v082.i13}
#'
#' @return An object of class `mlm_probe` (a list) with components:
#'   * `slopes`: a data frame with columns `modx_value`, `slope`, `se`, `t`,
#'     `df`, `p`, `ci_lower`, `ci_upper`.
#'   * `df_method`: the degrees-of-freedom method used.
#'   * `pred`, `modx`: names of the predictor and moderator.
#'   * `modx.values`: the strategy used.
#'   * `conf.level`: the confidence level.
#'   * `model`: the original model (stored for downstream use).
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200), x = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' mlm_probe(mod, pred = "x", modx = "m")
#' mlm_probe(mod, pred = "x", modx = "m", modx.values = "quartiles")
#' mlm_probe(mod, pred = "x", modx = "m", modx.values = "custom", at = c(-1, 0, 1))
#'
#' @export
mlm_probe <- function(model,
                      pred,
                      modx,
                      modx.values = c("mean-sd", "quartiles", "tertiles", "custom"),
                      at          = NULL,
                      conf.level  = 0.95,
                      df_method   = c("satterthwaite", "kenward-roger",
                                      "between", "residual"),
                      modx.level  = c("auto", "cluster", "observation")) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  modx.values <- match.arg(modx.values)
  df_method   <- match.arg(df_method)
  modx.level  <- match.arg(modx.level)

  ctx       <- .infer_ctx(model, df_method)
  modx_vals <- .pick_modx_values(.modx_vector(model, modx, modx.level),
                                 modx.values = modx.values, at = at)

  rows <- lapply(modx_vals, function(w) {
    ss <- .simple_slope(ctx, pred, modx, w, conf.level = conf.level)
    data.frame(modx_value = w, slope = ss$slope, se = ss$se, t = ss$t,
               df = ss$df, p = ss$p, ci_lower = ss$ci_lower,
               ci_upper = ss$ci_upper, stringsAsFactors = FALSE)
  })
  slopes_df <- do.call(rbind, rows)
  rownames(slopes_df) <- NULL

  structure(
    list(
      slopes      = slopes_df,
      pred        = pred,
      modx        = modx,
      modx.values = modx.values,
      conf.level  = conf.level,
      df_method   = df_method,
      model       = model
    ),
    class = "mlm_probe"
  )
}

#' @export
print.mlm_probe <- function(x, digits = 3, ...) {
  cat("\n--- Simple Slopes: mlm_probe ---\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Confidence level:", x$conf.level, "\n")
  if (!is.null(x$df_method)) cat("df method       :", x$df_method, "\n")
  cat("\n")

  df <- x$slopes
  df_print <- data.frame(
    modx_value = round(df$modx_value, digits),
    slope      = round(df$slope,      digits),
    se         = round(df$se,         digits),
    t          = round(df$t,          digits),
    df         = round(df$df,         1),
    p          = format_pval(df$p),
    ci_lower   = round(df$ci_lower,   digits),
    ci_upper   = round(df$ci_upper,   digits)
  )
  names(df_print)[1] <- x$modx

  print(df_print, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' Format p-values for printing
#' @noRd
format_pval <- function(p) {
  ifelse(p < 0.001, "< .001", formatC(p, digits = 3, format = "f"))
}
