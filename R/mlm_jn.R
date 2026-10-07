#' Johnson---Neyman interval for multilevel two-way interactions
#'
#' Computes the Johnson---Neyman (JN) interval: the region(s) of the moderator
#' (`modx`) where the simple slope of `pred` transitions between statistical
#' significance and non-significance. Useful for identifying *exactly* at which
#' moderator values an effect becomes significant.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx` in the fixed-effects structure.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param alpha Significance level. Default `0.05`.
#' @param modx.range Numeric vector of length 2 giving the range over which to
#'   evaluate the slope. Defaults to the observed range of `modx`.
#' @param grid Integer. Number of points at which the slope is evaluated
#'   across the moderator range (used for plotting and to bracket the
#'   boundaries). Default `200`.
#' @inheritParams mlm_probe
#'
#' @details
#' A Johnson-Neyman boundary is a moderator value \eqn{w} at which
#' \eqn{|t(w)| = t_{1-\alpha/2,\,df(w)}}, where
#' \eqn{t(w) = (\hat\beta_1 + \hat\beta_3 w) / SE(w)}. With a constant
#' number of degrees of freedom (`df_method = "between"` or `"residual"`) the
#' boundaries are the real roots of a quadratic in \eqn{w} (Bauer and Curran,
#' 2005), and these closed-form roots are returned. With Satterthwaite or
#' Kenward-Roger degrees of freedom, which change with \eqn{w}, the
#' boundaries are found by root-finding (`stats::uniroot()`) on the exact
#' criterion, starting from brackets located on the grid.
#'
#' Only boundaries inside `modx.range` are returned in `jn_bounds`. The
#' closed-form roots outside that range are reported in `jn_bounds_all`
#' when available, because a boundary outside the observed data is an
#' extrapolation.
#'
#' @return An object of class `mlm_jn` with components:
#'   * `jn_bounds`: numeric vector of moderator values inside `modx.range`
#'     where the slope crosses the significance threshold. `NA` if none.
#'   * `jn_bounds_all`: for constant-df methods, all real roots of the
#'     Johnson-Neyman quadratic (possibly outside the data); `NULL` otherwise.
#'   * `slopes_df`: data frame of slope estimates and significance across the
#'     grid.
#'   * `pred`, `modx`, `alpha`, `modx.range`, `df_method`.
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
#' jn <- mlm_jn(mod, pred = "x", modx = "m")
#' print(jn)
#'
#' @export
mlm_jn <- function(model,
                   pred,
                   modx,
                   alpha       = 0.05,
                   modx.range  = NULL,
                   grid        = 200L,
                   df_method   = c("satterthwaite", "kenward-roger",
                                   "between", "residual")) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  df_method <- match.arg(df_method)
  ctx <- .infer_ctx(model, df_method)

  modx_x <- model@frame[[modx]]
  if (is.null(modx.range)) modx.range <- range(modx_x, na.rm = TRUE)
  modx_grid <- seq(modx.range[1], modx.range[2], length.out = grid)

  slopes_list <- lapply(modx_grid, function(w) {
    ss <- .simple_slope(ctx, pred, modx, w, conf.level = 1 - alpha)
    data.frame(modx_value = w, slope = ss$slope, se = ss$se, t = ss$t,
               df = ss$df, p = ss$p, sig = ss$p < alpha,
               ci_lower = ss$ci_lower, ci_upper = ss$ci_upper,
               stringsAsFactors = FALSE)
  })
  slopes_df <- do.call(rbind, slopes_list)
  rownames(slopes_df) <- NULL

  # Criterion: |t(w)| - t_crit(df(w)); zero at a JN boundary.
  crit <- function(w) {
    ss <- .simple_slope(ctx, pred, modx, w, conf.level = 1 - alpha)
    abs(ss$t) - stats::qt(1 - alpha / 2, ss$df)
  }

  jn_all <- NULL
  if (df_method %in% c("between", "residual")) {
    jn_all <- .jn_closed_form(ctx, pred, modx, alpha)
    inside <- jn_all[jn_all >= modx.range[1] & jn_all <= modx.range[2]]
    jn_bounds <- if (length(inside)) sort(inside) else NA_real_
  } else {
    g  <- abs(slopes_df$t) - stats::qt(1 - alpha / 2, slopes_df$df)
    ix <- which(diff(sign(g)) != 0)
    jn_bounds <- if (length(ix) == 0L) NA_real_ else vapply(ix, function(i) {
      stats::uniroot(crit, c(modx_grid[i], modx_grid[i + 1]),
                     tol = 1e-10)$root
    }, numeric(1))
  }

  structure(
    list(
      jn_bounds     = jn_bounds,
      jn_bounds_all = jn_all,
      slopes_df     = slopes_df,
      pred          = pred,
      modx          = modx,
      alpha         = alpha,
      modx.range    = modx.range,
      df_method     = df_method
    ),
    class = "mlm_jn"
  )
}

# Closed-form JN roots (Bauer and Curran, 2005) for constant df.
.jn_closed_form <- function(ctx, pred, modx, alpha) {
  int_term <- .get_interaction_term(ctx$model, pred, modx)
  b1 <- ctx$fe[[pred]]; b3 <- ctx$fe[[int_term]]
  v11 <- ctx$V[pred, pred]; v33 <- ctx$V[int_term, int_term]
  v13 <- ctx$V[pred, int_term]
  tc2 <- stats::qt(1 - alpha / 2, ctx$df_const)^2
  a <- b3^2 - tc2 * v33
  b <- 2 * (b1 * b3 - tc2 * v13)
  c <- b1^2 - tc2 * v11
  if (abs(a) < .Machine$double.eps) {
    return(if (abs(b) < .Machine$double.eps) numeric(0) else -c / b)
  }
  disc <- b^2 - 4 * a * c
  if (disc < 0) return(numeric(0))
  sort(c((-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a)))
}

#' Plot a Johnson-Neyman interval from a multilevel model
#'
#' Creates a ggplot2 figure showing the simple slope of `pred` across the
#' full range of `modx`, with shading indicating regions of significance.
#' A vertical dashed line marks each Johnson-Neyman boundary.
#'
#' @param x An `mlm_jn` object from `mlm_jn()`.
#' @param x_label Label for the x-axis. Defaults to the moderator name.
#' @param y_label Label for the y-axis. Defaults to "Simple slope of pred".
#' @param sig_color Fill colour for the significant region. Default `"#2166AC"`.
#' @param nonsig_color Fill colour for the non-significant region. Default `"#D6604D"`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
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
#' jn <- mlm_jn(mod, pred = "x", modx = "m")
#' plot(jn)
#'
#' @export
plot.mlm_jn <- function(x,
                        x_label    = NULL,
                        y_label    = NULL,
                        sig_color  = "#2166AC",
                        nonsig_color = "#D6604D",
                        ...) {
  df <- x$slopes_df

  if (is.null(x_label)) x_label <- x$modx
  if (is.null(y_label)) y_label <- paste0("Simple slope of '", x$pred, "'")

  # Significance as factor for fill
  df$region <- ifelse(df$sig, "Significant", "Non-significant")
  df$region <- factor(df$region, levels = c("Significant", "Non-significant"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = modx_value, y = slope)) +
    # Shaded ribbon for CI
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper, fill = region),
      alpha = 0.25
    ) +
    # Slope line coloured by significance
    ggplot2::geom_line(ggplot2::aes(color = region), linewidth = 1.1) +
    # Zero reference line
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.7) +
    # JN boundary lines
    {
      if (!all(is.na(x$jn_bounds))) {
        lapply(x$jn_bounds, function(b) {
          ggplot2::geom_vline(xintercept = b, linetype = "dashed",
                              color = "grey20", linewidth = 0.8)
        })
      }
    } +
    # JN boundary labels
    {
      if (!all(is.na(x$jn_bounds))) {
        lapply(x$jn_bounds, function(b) {
          ggplot2::annotate("text",
            x = b, y = max(df$ci_upper, na.rm = TRUE),
            label = paste0(x$modx, " = ", round(b, 2)),
            hjust = -0.05, vjust = 1, size = 3.5, color = "grey20"
          )
        })
      }
    } +
    ggplot2::scale_color_manual(
      values = c("Significant" = sig_color, "Non-significant" = nonsig_color),
      name   = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Significant" = sig_color, "Non-significant" = nonsig_color),
      name   = NULL
    ) +
    ggplot2::labs(x = x_label, y = y_label,
                  title = "Johnson-Neyman Plot",
                  subtitle = paste0("Regions where slope of '", x$pred,
                                    "' is significant (p < ", x$alpha, ")")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = 11)
    )

  p
}

#' @export
print.mlm_jn <- function(x, digits = 3, ...) {
  cat("\n--- Johnson-Neyman Interval: mlm_jn ---\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Alpha           :", x$alpha, "\n")
  cat("Moderator range :", round(x$modx.range[1], digits), "to",
      round(x$modx.range[2], digits), "\n\n")

  if (all(is.na(x$jn_bounds))) {
    cat("No Johnson-Neyman boundary found in the observed moderator range.\n")
    # Indicate whether slope is universally sig or non-sig
    pct_sig <- mean(x$slopes_df$sig, na.rm = TRUE)
    if (pct_sig > 0.99) {
      cat("The simple slope of", x$pred, "is significant across the entire",
          "moderator range.\n")
    } else if (pct_sig < 0.01) {
      cat("The simple slope of", x$pred, "is non-significant across the entire",
          "moderator range.\n")
    }
  } else {
    cat("Johnson-Neyman boundary/boundaries (", x$modx, "):\n")
    for (b in x$jn_bounds) {
      cat("  ", round(b, digits), "\n")
    }
    cat("\nInterpretation: The simple slope of '", x$pred, "' is statistically\n",
        "significant (p < ", x$alpha, ") for values of '", x$modx, "'\n", sep = "")

    # Describe each significant stretch between consecutive boundaries.
    cuts <- c(x$modx.range[1], sort(x$jn_bounds), x$modx.range[2])
    segs <- character(0)
    for (k in seq_len(length(cuts) - 1L)) {
      mid <- (cuts[k] + cuts[k + 1]) / 2
      i <- which.min(abs(x$slopes_df$modx_value - mid))
      if (isTRUE(x$slopes_df$sig[i])) {
        segs <- c(segs, sprintf("%s to %s", round(cuts[k], digits),
                                round(cuts[k + 1], digits)))
      }
    }
    if (length(segs)) {
      cat("  in ", paste(segs, collapse = " and "),
          " (within the observed range).\n", sep = "")
    } else {
      cat("  nowhere in the observed range.\n")
    }
  }
  if (!is.null(x$df_method)) cat("df method:", x$df_method, "\n")
  cat("\n")
  invisible(x)
}
