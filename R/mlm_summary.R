#' Summary table for a multilevel moderation effect
#'
#' Returns a consolidated summary of the moderation effect: the focal
#' interaction coefficient, simple slopes at selected moderator values, and
#' (optionally) the Johnson---Neyman interval. Designed for quick reporting and
#' results sections.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param modx.values Moderator value strategy. See `mlm_probe()`.
#' @param at Optional numeric vector of custom moderator values.
#' @param conf.level Confidence level. Default `0.95`.
#' @param jn Logical. Include Johnson-Neyman region? Default `TRUE`.
#' @param alpha Alpha for JN interval. Default `0.05`.
#' @inheritParams mlm_probe
#'
#' @return An object of class `mlm_summary` (a list) with components:
#'   * `interaction`: one-row data frame for the interaction term.
#'   * `simple_slopes`: data frame from `mlm_probe()`.
#'   * `jn`: output of `mlm_jn()` (or `NULL` if `jn = FALSE`).
#'   * Other metadata.
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
#' mlm_summary(mod, pred = "x", modx = "m")
#'
#' @export
mlm_summary <- function(model,
                        pred,
                        modx,
                        modx.values = c("mean-sd", "quartiles", "tertiles", "custom"),
                        at          = NULL,
                        conf.level  = 0.95,
                        jn          = TRUE,
                        alpha       = 0.05,
                        df_method   = c("satterthwaite", "kenward-roger",
                                        "between", "residual"),
                        modx.level  = c("auto", "cluster", "observation")) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  modx.values <- match.arg(modx.values)
  df_method   <- match.arg(df_method)
  modx.level  <- match.arg(modx.level)

  ctx      <- .infer_ctx(model, df_method)
  int_term <- .get_interaction_term(model, pred, modx)
  L <- stats::setNames(numeric(length(ctx$fe)), names(ctx$fe))
  L[int_term] <- 1
  r      <- .lincomb(ctx, L)
  t_int  <- r$estimate / r$se
  p_int  <- 2 * stats::pt(abs(t_int), df = r$df, lower.tail = FALSE)
  t_crit <- stats::qt(1 - (1 - conf.level) / 2, df = r$df)

  interaction_row <- data.frame(
    term     = int_term,
    estimate = r$estimate,
    se       = r$se,
    t        = t_int,
    df       = r$df,
    p        = p_int,
    ci_lower = r$estimate - t_crit * r$se,
    ci_upper = r$estimate + t_crit * r$se,
    stringsAsFactors = FALSE
  )

  probe_out <- mlm_probe(model, pred = pred, modx = modx,
                         modx.values = modx.values, at = at,
                         conf.level = conf.level, df_method = df_method,
                         modx.level = modx.level)

  jn_out <- NULL
  if (jn) {
    jn_out <- mlm_jn(model, pred = pred, modx = modx, alpha = alpha,
                     df_method = df_method)
  }

  structure(
    list(
      interaction   = interaction_row,
      simple_slopes = probe_out$slopes,
      jn            = jn_out,
      pred          = pred,
      modx          = modx,
      modx.values   = modx.values,
      conf.level    = conf.level,
      alpha         = alpha,
      df_method     = df_method,
      model         = model
    ),
    class = "mlm_summary"
  )
}

#' @export
print.mlm_summary <- function(x, digits = 3, ...) {
  cat("\n========================================\n")
  cat("  Multilevel Moderation Summary\n")
  cat("========================================\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Confidence level:", x$conf.level, "\n")
  if (!is.null(x$df_method)) cat("df method       :", x$df_method, "\n")

  cat("\n--- Interaction Term ---\n")
  ir <- x$interaction
  cat(sprintf(
    "  %-28s  b = %7.3f  SE = %6.3f  t(%g) = %6.3f  p = %s  [%s, %s]\n",
    ir$term,
    ir$estimate, ir$se, round(ir$df, 1), ir$t,
    format_pval(ir$p),
    round(ir$ci_lower, digits),
    round(ir$ci_upper, digits)
  ))

  cat("\n--- Simple Slopes ---\n")
  df <- x$simple_slopes
  header <- sprintf("  %-10s  %8s  %8s  %8s  %8s  %8s  %8s  %8s",
                    x$modx, "slope", "se", "t", "df", "p",
                    paste0(round(x$conf.level * 100), "% CI lo"),
                    paste0(round(x$conf.level * 100), "% CI hi"))
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf("  %-10.3f  %8.3f  %8.3f  %8.3f  %8.1f  %8s  %8.3f  %8.3f\n",
                df$modx_value[i], df$slope[i], df$se[i], df$t[i],
                df$df[i], format_pval(df$p[i]),
                df$ci_lower[i], df$ci_upper[i]))
  }

  if (!is.null(x$jn)) {
    cat("\n--- Johnson-Neyman Region ---\n")
    jb <- x$jn$jn_bounds
    if (all(is.na(jb))) {
      pct_sig <- mean(x$jn$slopes_df$sig, na.rm = TRUE)
      if (pct_sig > 0.99) {
        cat("  Slope of '", x$pred, "' is significant across the entire range of '",
            x$modx, "'.\n", sep = "")
      } else if (pct_sig < 0.01) {
        cat("  Slope of '", x$pred, "' is non-significant across the entire range of '",
            x$modx, "'.\n", sep = "")
      } else {
        cat("  No clear JN boundary found.\n")
      }
    } else {
      cat("  JN boundary/boundaries for '", x$modx, "': ",
          paste(round(jb, digits), collapse = ", "), "\n", sep = "")
    }
  }

  cat("\n========================================\n\n")
  invisible(x)
}
