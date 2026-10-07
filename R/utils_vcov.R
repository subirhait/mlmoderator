# ---------------------------------------------------------------------------
# Inference engine
#
# All tests and intervals in mlmoderator go through .lincomb(), which
# returns the estimate, standard error, and denominator degrees of freedom
# for a linear combination L'beta of the fixed effects. The degrees of
# freedom depend on `df_method`:
#
#   "satterthwaite" (default)  lmerTest::contest1D(), Satterthwaite (1946)
#   "kenward-roger"            lmerTest::contest1D() via pbkrtest; this also
#                              uses the Kenward-Roger adjusted covariance
#   "between"                  J - q - 1, where J is the number of clusters
#                              and q the number of cluster-level fixed
#                              effects (Snijders and Bosker, 2012)
#   "residual"                 N - p (the rule used before 0.3.0; it treats
#                              level-1 rows as independent and is
#                              anti-conservative for cross-level terms)
# ---------------------------------------------------------------------------

.df_methods <- c("satterthwaite", "kenward-roger", "between", "residual")

.extract_vcov <- function(model) {
  as.matrix(stats::vcov(model))
}

# Grouping factor used for cluster-level quantities: the first one.
.cluster_name <- function(model) {
  names(lme4::getME(model, "flist"))[1]
}

.cluster_id <- function(model) {
  lme4::getME(model, "flist")[[1]]
}

# Number of fixed-effect columns that are constant within every cluster
# (cluster-level predictors), excluding the intercept.
.n_level2_terms <- function(model) {
  X  <- lme4::getME(model, "X")
  id <- .cluster_id(model)
  keep <- colnames(X) != "(Intercept)"
  if (!any(keep)) return(0L)
  const <- vapply(which(keep), function(j) {
    all(tapply(X[, j], id, function(v) diff(range(v)) < 1e-12))
  }, logical(1))
  sum(const)
}

.df_between <- function(model) {
  J <- nlevels(.cluster_id(model))
  max(J - .n_level2_terms(model) - 1L, 1L)
}

.df_residual <- function(model) {
  max(nrow(model@frame) - length(lme4::fixef(model)), 1L)
}

# Kept for backward compatibility with code that called it directly.
get_residual_df <- function(model) .df_residual(model)

# Build an inference context once per call.
.infer_ctx <- function(model, df_method = "satterthwaite") {
  df_method <- match.arg(df_method, .df_methods)
  ctx <- list(model = model, method = df_method,
              fe = lme4::fixef(model), V = .extract_vcov(model))
  if (df_method %in% c("satterthwaite", "kenward-roger")) {
    if (df_method == "kenward-roger" &&
        !requireNamespace("pbkrtest", quietly = TRUE)) {
      rlang::abort("df_method = 'kenward-roger' requires the 'pbkrtest' package.")
    }
    ctx$mt <- tryCatch(
      suppressMessages(lmerTest::as_lmerModLmerTest(model)),
      error = function(e) rlang::abort(paste0(
        "Could not prepare the model for ", df_method,
        " degrees of freedom: ", conditionMessage(e),
        "\nTry df_method = 'between'.")))
  } else if (df_method == "between") {
    ctx$df_const <- .df_between(model)
  } else {
    ctx$df_const <- .df_residual(model)
  }
  ctx
}

# Estimate, SE, and df of L'beta.
.lincomb <- function(ctx, L) {
  L <- as.numeric(L)
  if (ctx$method %in% c("satterthwaite", "kenward-roger")) {
    ddf <- if (ctx$method == "satterthwaite") "Satterthwaite" else "Kenward-Roger"
    r <- lmerTest::contest1D(ctx$mt, L, ddf = ddf)
    return(list(estimate = unname(r[["Estimate"]]),
                se = unname(r[["Std. Error"]]),
                df = unname(r[["df"]])))
  }
  est <- sum(L * ctx$fe)
  se  <- sqrt(as.numeric(crossprod(L, ctx$V %*% L)))
  list(estimate = est, se = se, df = ctx$df_const)
}

# L vector for the simple slope of `pred` at moderator value w.
.slope_L <- function(ctx, pred, int_term, w) {
  L <- stats::setNames(numeric(length(ctx$fe)), names(ctx$fe))
  L[pred] <- 1
  L[int_term] <- w
  L
}

.simple_slope <- function(ctx, pred, modx, w, conf.level = 0.95) {
  int_term <- .get_interaction_term(ctx$model, pred, modx)
  r <- .lincomb(ctx, .slope_L(ctx, pred, int_term, w))
  t_val <- r$estimate / r$se
  tc <- stats::qt(1 - (1 - conf.level) / 2, df = r$df)
  list(slope = r$estimate, se = r$se, t = t_val, df = r$df,
       p = 2 * stats::pt(abs(t_val), df = r$df, lower.tail = FALSE),
       ci_lower = r$estimate - tc * r$se,
       ci_upper = r$estimate + tc * r$se)
}

# Backward-compatible wrapper (0.2.x internal name).
.simple_slope_linear <- function(model, pred, modx, modx_val,
                                 conf.level = 0.95, df_method = "satterthwaite") {
  .simple_slope(.infer_ctx(model, df_method), pred, modx, modx_val, conf.level)
}

# Fixed-effects design matrix for new data (for fitted-value SEs).
.fixed_X <- function(model, newdata) {
  tt <- stats::delete.response(stats::terms(model, fixed.only = TRUE))
  xl <- .xlevels(model)
  xl <- xl[names(xl) %in% all.vars(tt)]
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                           xlev = if (length(xl)) xl else NULL)
  X  <- stats::model.matrix(tt, mf,
                            contrasts.arg = attr(lme4::getME(model, "X"), "contrasts"))
  X[, names(lme4::fixef(model)), drop = FALSE]
}

.xlevels <- function(model) {
  fr <- model@frame
  f  <- vapply(fr, is.factor, logical(1))
  lapply(fr[f], levels)
}
