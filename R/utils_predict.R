# Internal prediction helpers

# Moderator values used to choose probe points. For a moderator that is
# constant within clusters (a cluster-level moderator), the default uses one
# value per cluster so that the mean and SD describe clusters rather than
# being weighted by cluster size.
.modx_vector <- function(model, modx, modx.level = c("auto", "cluster", "observation")) {
  modx.level <- match.arg(modx.level)
  x  <- model@frame[[modx]]
  id <- .cluster_id(model)
  is_l2 <- all(tapply(x, id, function(v) diff(range(v)) < 1e-12))
  if (modx.level == "observation" || (modx.level == "auto" && !is_l2)) return(x)
  if (modx.level == "cluster" && !is_l2) {
    rlang::abort(paste0("'", modx, "' varies within clusters; ",
                        "use modx.level = 'observation'."))
  }
  as.numeric(tapply(x, id, function(v) v[1]))
}

.pick_modx_values <- function(x, modx.values = "mean-sd", at = NULL) {
  # Custom values via `at` always win, regardless of modx.values
  if (!is.null(at)) return(sort(as.numeric(at)))

  modx.values <- match.arg(
    modx.values,
    choices = c("mean-sd", "quartiles", "tertiles", "custom")
  )

  switch(modx.values,
    "mean-sd" = {
      m <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      c(m - s, m, m + s)
    },
    "quartiles" = {
      as.numeric(stats::quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
    },
    "tertiles" = {
      as.numeric(stats::quantile(x, probs = c(1/3, 2/3), na.rm = TRUE))
    },
    "custom" = {
      rlang::abort(
        "Supply moderator values via the `at` argument when modx.values = 'custom'."
      )
    }
  )
}

.make_prediction_grid <- function(model, pred, modx, modx_vals, n_pred = 100) {
  mf <- model@frame
  cluster_vars <- names(lme4::getME(model, "flist"))

  pred_range <- seq(
    min(mf[[pred]], na.rm = TRUE),
    max(mf[[pred]], na.rm = TRUE),
    length.out = n_pred
  )

  grids <- lapply(modx_vals, function(w) {
    g <- mf[rep(1L, n_pred), , drop = FALSE]
    rownames(g) <- NULL
    g[[pred]] <- pred_range
    g[[modx]] <- w
    g <- .hold_covariates(g, mf, c(pred, modx, cluster_vars))
    g$.modx_val <- w
    g
  })

  do.call(rbind, grids)
}

# Hold every other variable at its mean (numeric) or reference level (factor).
.hold_covariates <- function(g, mf, exclude) {
  other_vars <- setdiff(names(g), exclude)
  for (v in other_vars) {
    if (is.numeric(g[[v]])) {
      g[[v]] <- mean(mf[[v]], na.rm = TRUE)
    } else if (is.factor(g[[v]])) {
      g[[v]] <- factor(levels(mf[[v]])[1], levels = levels(mf[[v]]))
    } else if (is.character(g[[v]])) {
      g[[v]] <- mf[[v]][1]
    }
  }
  g
}

.build_modx_labels <- function(vals, modx, strategy) {
  lbls <- if (strategy == "mean-sd" && length(vals) == 3L) {
    c("-1 SD", "Mean", "+1 SD")
  } else if (strategy == "quartiles") {
    c("25th pct", "50th pct", "75th pct")
  } else if (strategy == "tertiles") {
    c("1st tertile", "2nd tertile")
  } else {
    as.character(round(vals, 2))
  }
  stats::setNames(lbls, as.character(vals))
}
