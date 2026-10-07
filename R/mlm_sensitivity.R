#' Leave-one-cluster-out influence diagnostics for a cross-level interaction
#'
#' Refits the model once for each cluster, leaving that cluster out, and
#' records how the interaction coefficient changes. This answers a specific
#' question: *is the interaction driven by a small number of clusters?*
#'
#' Each refit uses [stats::update()] on the data the model was fitted to
#' (via [lme4::getData()]), so it keeps the original formula, estimation
#' method (REML or ML), weights, and control settings. Influence is measured
#' by a DFBETA-type statistic,
#' \deqn{\mathrm{DFBETA}_j = (\hat\beta_{3(-j)} - \hat\beta_3) / SE(\hat\beta_3),}
#' the change in the interaction when cluster \eqn{j} is omitted, in units of
#' the full-data standard error. Clusters with
#' \eqn{|\mathrm{DFBETA}_j| > 2/\sqrt{J}} are flagged, following the usual
#' size-adjusted cutoff for DFBETAS (Belsley, Kuh, and Welsch, 1980). The
#' cutoff is a screening heuristic, not a test.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param alpha Significance level used to record whether each refit keeps
#'   the interaction significant. Default `0.05`.
#' @param verbose Logical. Print progress during refitting? Default `FALSE`.
#' @param icc_range,icc_grid,loco Defunct arguments from versions before
#'   0.3.0 (see Details). Supplying `icc_range` or `icc_grid` gives a warning
#'   and is otherwise ignored; `loco = FALSE` is no longer supported.
#' @inheritParams mlm_probe
#'
#' @details
#' **Changes in 0.3.0.** Earlier versions also reported an "ICC-shift"
#' analysis that rescaled the interaction's standard error by a
#' design-effect ratio for hypothetical intraclass correlations. That formula
#' applies to means under a random-intercept model, not to a cross-level
#' interaction in a random-slope model, and the accompanying adjusted
#' Johnson-Neyman boundary was not correct, so the analysis has been removed.
#' The earlier leave-one-cluster-out refits used maximum likelihood even when
#' the original model used REML, and refitted from the model frame, which
#' failed silently for formulas with transformed variables; both are fixed.
#'
#' **Scope.** These are influence diagnostics. They do not address
#' unmeasured confounding of the interaction.
#'
#' @references
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). *Regression diagnostics*.
#' Wiley.
#'
#' @return An object of class `mlm_sensitivity` with components:
#'   * `observed`: list with the full-data interaction estimate, SE, df,
#'     p-value, and number of clusters.
#'   * `loco`: data frame with one row per cluster: `cluster`, `n_obs`,
#'     `b_int`, `se_int`, `df`, `p_int`, `b_change`, `dfbeta`, `influential`,
#'     `same_sign`, `sig`, and `error` (the error message when a refit
#'     failed, otherwise `NA`).
#'   * `cutoff`: the DFBETA flagging cutoff \eqn{2/\sqrt{J}}.
#'   * `summary`: list with the number of failed refits, the proportion of
#'     successful refits keeping the sign and significance of the interaction,
#'     and the range of refitted estimates.
#'   * Metadata: `pred`, `modx`, `alpha`, `df_method`, `grp_name`,
#'     `int_term`.
#'
#' @examples
#' set.seed(42)
#' n_j <- 20; n_i <- 10
#' dat_small <- data.frame(
#'   y   = rnorm(n_j * n_i),
#'   x   = rnorm(n_j * n_i),
#'   m   = rep(rnorm(n_j), each = n_i),
#'   grp = factor(rep(seq_len(n_j), each = n_i))
#' )
#' dat_small$y <- dat_small$y + 0.5 * dat_small$x * dat_small$m
#' mod_small <- lme4::lmer(y ~ x * m + (1 | grp), data = dat_small)
#'
#' sens <- mlm_sensitivity(mod_small, pred = "x", modx = "m",
#'                         df_method = "between")
#' sens
#' plot(sens)
#'
#' @export
mlm_sensitivity <- function(model,
                            pred,
                            modx,
                            alpha     = 0.05,
                            df_method = c("satterthwaite", "kenward-roger",
                                          "between", "residual"),
                            verbose   = FALSE,
                            icc_range,
                            icc_grid,
                            loco      = TRUE) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  df_method <- match.arg(df_method)

  if (!missing(icc_range) || !missing(icc_grid)) {
    rlang::warn(paste0(
      "The ICC-shift analysis was removed in mlmoderator 0.3.0; ",
      "`icc_range` and `icc_grid` are ignored. See ?mlm_sensitivity."))
  }
  if (!isTRUE(loco)) {
    rlang::abort(paste0(
      "`loco = FALSE` is no longer supported: leave-one-cluster-out ",
      "refitting is the only analysis mlm_sensitivity() performs."))
  }

  int_term <- .get_interaction_term(model, pred, modx)
  grp_name <- .cluster_name(model)

  ctx0 <- .infer_ctx(model, df_method)
  L0   <- stats::setNames(numeric(length(ctx0$fe)), names(ctx0$fe))
  L0[int_term] <- 1
  r0 <- .lincomb(ctx0, L0)
  p0 <- 2 * stats::pt(abs(r0$estimate / r0$se), r0$df, lower.tail = FALSE)

  dat <- tryCatch(lme4::getData(model), error = function(e) NULL)
  if (is.null(dat) || !is.data.frame(dat)) {
    rlang::abort(paste0(
      "Could not recover the data the model was fitted to. Fit the model ",
      "with a `data =` argument that is available in the calling environment."))
  }
  if (!grp_name %in% names(dat)) {
    rlang::abort(paste0(
      "Grouping factor '", grp_name, "' is not a column of the model data. ",
      "Nested or constructed grouping factors (e.g. a:b) are not supported; ",
      "create the grouping variable as a column first."))
  }

  cl_ids <- unique(as.character(dat[[grp_name]]))
  cl_ids <- cl_ids[!is.na(cl_ids)]
  J <- length(cl_ids)
  if (verbose) cat("Refitting without each of", J, "clusters...\n")

  rows <- lapply(seq_along(cl_ids), function(i) {
    cid <- cl_ids[i]
    if (verbose && i %% 10 == 0) cat("  cluster", i, "of", J, "\n")
    keep <- as.character(dat[[grp_name]]) != cid | is.na(dat[[grp_name]])
    dsub <- dat[keep, , drop = FALSE]
    err  <- NA_character_
    fit  <- tryCatch(
      suppressMessages(suppressWarnings(stats::update(model, data = dsub))),
      error = function(e) { err <<- conditionMessage(e); NULL })
    if (!is.null(fit) && !int_term %in% names(lme4::fixef(fit))) {
      err <- "interaction not estimable"; fit <- NULL
    }
    if (is.null(fit)) {
      return(data.frame(cluster = cid, n_obs = sum(!keep), b_int = NA_real_,
                        se_int = NA_real_, df = NA_real_, p_int = NA_real_,
                        error = err, stringsAsFactors = FALSE))
    }
    ctx <- tryCatch(.infer_ctx(fit, df_method),
                    error = function(e) .infer_ctx(fit, "between"))
    L <- stats::setNames(numeric(length(ctx$fe)), names(ctx$fe))
    L[int_term] <- 1
    r <- .lincomb(ctx, L)
    data.frame(cluster = cid, n_obs = sum(!keep), b_int = r$estimate,
               se_int = r$se, df = r$df,
               p_int = 2 * stats::pt(abs(r$estimate / r$se), r$df,
                                     lower.tail = FALSE),
               error = NA_character_, stringsAsFactors = FALSE)
  })
  loco_df <- do.call(rbind, rows)
  rownames(loco_df) <- NULL

  cutoff <- 2 / sqrt(J)
  loco_df$b_change    <- loco_df$b_int - r0$estimate
  loco_df$dfbeta      <- loco_df$b_change / r0$se
  loco_df$influential <- abs(loco_df$dfbeta) > cutoff
  loco_df$same_sign   <- sign(loco_df$b_int) == sign(r0$estimate)
  loco_df$sig         <- loco_df$p_int < alpha
  loco_df <- loco_df[, c("cluster", "n_obs", "b_int", "se_int", "df", "p_int",
                         "b_change", "dfbeta", "influential", "same_sign",
                         "sig", "error")]

  n_fail <- sum(is.na(loco_df$b_int))
  if (n_fail > 0) {
    rlang::warn(sprintf(
      "%d of %d leave-one-cluster-out refits failed; see the `error` column.",
      n_fail, J))
  }
  ok <- !is.na(loco_df$b_int)
  keeps <- if (any(ok)) {
    mean(loco_df$same_sign[ok] & (loco_df$sig[ok] == (p0 < alpha)))
  } else NA_real_

  structure(
    list(
      observed = list(b_int = r0$estimate, se_int = r0$se, df = r0$df,
                      p_int = p0, n_clusters = J, n_total = nrow(model@frame)),
      loco     = loco_df,
      cutoff   = cutoff,
      summary  = list(
        n_failed          = n_fail,
        prop_same_result  = keeps,
        b_range           = if (any(ok)) range(loco_df$b_int[ok]) else c(NA, NA),
        n_influential     = sum(loco_df$influential, na.rm = TRUE)
      ),
      pred      = pred,
      modx      = modx,
      alpha     = alpha,
      df_method = df_method,
      grp_name  = grp_name,
      int_term  = int_term
    ),
    class = "mlm_sensitivity"
  )
}

#' @export
print.mlm_sensitivity <- function(x, digits = 3, ...) {
  o <- x$observed
  s <- x$summary
  cat("\n--- Leave-one-cluster-out influence: mlm_sensitivity ---\n")
  cat("Interaction     :", x$int_term, "\n")
  cat("Clusters        :", o$n_clusters, "(", x$grp_name, ")\n")
  cat("df method       :", x$df_method, "\n\n")
  cat(sprintf("Full data       : b = %.*f, SE = %.*f, df = %.1f, p = %s\n",
              digits, o$b_int, digits, o$se_int, o$df, format_pval(o$p_int)))
  cat(sprintf("Refitted range  : %.*f to %.*f\n", digits, s$b_range[1],
              digits, s$b_range[2]))
  cat(sprintf("Same sign and significance decision in %.1f%% of refits\n",
              100 * s$prop_same_result))
  cat(sprintf("Flagged (|DFBETA| > %.3f): %d cluster(s)\n", x$cutoff,
              s$n_influential))
  if (s$n_influential > 0) {
    fl <- x$loco[which(x$loco$influential), ]
    fl <- fl[order(-abs(fl$dfbeta)), ]
    for (i in seq_len(min(5L, nrow(fl)))) {
      cat(sprintf("  %-12s DFBETA = %6.3f  b without it = %.*f\n",
                  fl$cluster[i], fl$dfbeta[i], digits, fl$b_int[i]))
    }
  }
  if (s$n_failed > 0) cat("Failed refits   :", s$n_failed, "\n")
  cat("\n")
  invisible(x)
}

#' Plot leave-one-cluster-out influence for a cross-level interaction
#'
#' Plots the DFBETA of each cluster, ordered by size, with the flagging
#' cutoff \eqn{\pm 2/\sqrt{J}} as dashed lines.
#'
#' @param x An `mlm_sensitivity` object.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.mlm_sensitivity <- function(x, ...) {
  d <- x$loco[!is.na(x$loco$dfbeta), ]
  d <- d[order(d$dfbeta), ]
  d$cluster <- factor(d$cluster, levels = d$cluster)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$cluster, y = .data$dfbeta,
                                  color = .data$influential)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50") +
    ggplot2::geom_hline(yintercept = c(-1, 1) * x$cutoff,
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey30",
                                           `TRUE` = "#D6604D"),
                                name = "Flagged") +
    ggplot2::labs(
      x = paste0("Cluster omitted (", x$grp_name, "), ordered"),
      y = "DFBETA of the interaction",
      title = paste0("Leave-one-cluster-out influence: ", x$int_term)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}
