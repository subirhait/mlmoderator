#' Contour plot of predicted outcomes over the predictor x moderator space
#'
#' Plots iso-outcome contour lines of \eqn{\hat{Y}(x, w)} over the full joint
#' space of `pred` (x-axis) and `modx` (y-axis). This is the most direct
#' geometric representation of a two-way interaction:
#'
#' * **No interaction**: contour lines are perfectly straight and parallel ---
#'   the effect of `pred` does not depend on `modx`.
#' * **Positive interaction**: contour lines fan outward (rotate clockwise) ---
#'   higher `modx` steepens the `pred` slope.
#' * **Negative interaction**: contour lines fan inward (rotate counter-clockwise).
#'
#' The degree of non-parallelism among contours is a direct visual index of
#' interaction strength: the larger \eqn{\beta_3}, the more the lines rotate.
#'
#' An optional overlay draws the three standard simple-slope evaluation lines
#' (mean -- 1 SD of `modx`) as horizontal reference lines, connecting the plot
#' to `mlm_probe()` output.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx` in the fixed effects.
#' @param pred Character scalar. Focal predictor (x-axis).
#' @param modx Character scalar. Moderator (y-axis).
#' @param grid Integer. Grid resolution (points per axis). Default `80`.
#' @param n_contours Integer. Number of contour levels to draw. Default `10`.
#' @param fill Logical. Fill contour bands with colour? Default `TRUE`.
#' @param probe_lines Logical. Overlay horizontal lines at mean -- 1 SD of
#'   `modx`? Default `TRUE`.
#' @param x_label x-axis label. Defaults to `pred`.
#' @param y_label y-axis label. Defaults to `modx`.
#' @param legend_title Legend title. Defaults to the outcome variable name.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Predicted values are computed from fixed effects only
#' (`re.form = NA`), with all covariates held at their means or reference
#' levels. The surface therefore represents the population-average predicted
#' outcome, not any specific cluster.
#'
#' **Reading the plot:**
#' Pick any contour line. Its slope in the `pred` direction tells you how fast
#' the outcome changes with `pred` at that `modx` value. If the contour slopes
#' steeply up-right, `pred` has a strong positive effect there. If contours
#' become more horizontal as `modx` increases, the `pred` effect is weakening.
#' If they rotate from positive to flat to negative, you have a sign-changing
#' interaction --- and the `modx` value where they are perfectly horizontal is
#' the Johnson-Neyman boundary.
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
#' mlm_surface(mod, pred = "x", modx = "m")
#' mlm_surface(mod, pred = "x", modx = "m", fill = FALSE, n_contours = 15)
#'
#' @export
mlm_surface <- function(model,
                        pred,
                        modx,
                        grid         = 80L,
                        n_contours   = 10L,
                        fill         = TRUE,
                        probe_lines  = TRUE,
                        x_label      = NULL,
                        y_label      = NULL,
                        legend_title = NULL) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)

  mf      <- model@frame
  outcome <- names(mf)[1]

  if (is.null(x_label))      x_label      <- pred
  if (is.null(y_label))      y_label      <- modx
  if (is.null(legend_title)) legend_title <- outcome

  # ------ Build prediction grid ------------------------------------------------------------------------------------------------------------------------------------------------------
  pred_seq <- seq(min(mf[[pred]], na.rm = TRUE),
                  max(mf[[pred]], na.rm = TRUE),
                  length.out = grid)
  modx_seq <- seq(min(mf[[modx]], na.rm = TRUE),
                  max(mf[[modx]], na.rm = TRUE),
                  length.out = grid)

  grid_df <- expand.grid(pred_val = pred_seq, modx_val = modx_seq)

  # Build a newdata frame with all model variables
  template  <- mf[rep(1L, nrow(grid_df)), , drop = FALSE]
  rownames(template) <- NULL

  template[[pred]] <- grid_df$pred_val
  template[[modx]] <- grid_df$modx_val

  # Hold covariates at mean / reference level
  cluster_vars <- names(lme4::getME(model, "flist"))
  other_vars   <- setdiff(names(template), c(pred, modx, cluster_vars))
  for (v in other_vars) {
    if (is.numeric(template[[v]])) {
      template[[v]] <- mean(mf[[v]], na.rm = TRUE)
    } else if (is.factor(template[[v]])) {
      template[[v]] <- factor(levels(mf[[v]])[1], levels = levels(mf[[v]]))
    }
  }

  grid_df$yhat <- stats::predict(model, newdata = template, re.form = NA)

  # ------ Moderator probe values (mean -- 1 SD) ------------------------------------------------------------------------------------------------------
  modx_v <- .modx_vector(model, modx, "auto")
  m_modx <- mean(modx_v, na.rm = TRUE)
  s_modx <- stats::sd(modx_v, na.rm = TRUE)
  probe_vals <- c(m_modx - s_modx, m_modx, m_modx + s_modx)
  probe_labels <- c("-1 SD", "Mean", "+1 SD")

  # ------ Contour breaks ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  yhat_range <- range(grid_df$yhat, na.rm = TRUE)
  breaks <- pretty(yhat_range, n = n_contours)

  # ------ Plot ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  p <- ggplot2::ggplot(
    grid_df,
    ggplot2::aes(x = pred_val, y = modx_val, z = yhat)
  )

  if (fill) {
    p <- p +
      ggplot2::geom_contour_filled(
        breaks = breaks,
        alpha  = 0.85
      ) +
      ggplot2::scale_fill_viridis_d(
        name   = legend_title,
        option = "D",
        direction = 1
      )
  }

  p <- p +
    ggplot2::geom_contour(
      breaks    = breaks,
      color     = if (fill) "white" else "grey30",
      linewidth = if (fill) 0.3 else 0.6,
      alpha     = if (fill) 0.7 else 1
    )

  # Label contour lines
  p <- p +
    ggplot2::geom_contour(
      breaks    = breaks,
      color     = "white",
      linewidth = 0,
      alpha     = 0
    )

  # Probe lines (horizontal dashes at mean -- 1 SD of modx)
  if (probe_lines) {
    probe_df <- data.frame(
      y     = probe_vals,
      label = probe_labels
    )
    p <- p +
      ggplot2::geom_hline(
        data     = probe_df,
        ggplot2::aes(yintercept = y),
        linetype = "dashed",
        color    = "white",
        linewidth = 0.7,
        alpha    = 0.9
      ) +
      ggplot2::annotate(
        "text",
        x     = max(pred_seq),
        y     = probe_vals,
        label = probe_labels,
        hjust = 1.05, vjust = -0.5,
        size  = 3.2,
        color = "white",
        fontface = "bold"
      )
  }

  # Interaction annotation
  fe       <- lme4::fixef(model)
  int_term <- .get_interaction_term(model, pred, modx)
  b3       <- fe[int_term]
  b1       <- fe[pred]

  subtitle_text <- paste0(
    "dY/d", pred, " = ",
    round(b1, 3), " + ", round(b3, 3), " * ", modx,
    "   |   parallel contours = no interaction; rotation = moderation"
  )

  p +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = paste0("Predicted ", outcome,
                        ": contour surface over ", pred, " x ", modx),
      subtitle = subtitle_text
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.key.height = ggplot2::unit(1.4, "cm"),
      plot.subtitle     = ggplot2::element_text(color = "grey40", size = 9),
      panel.background  = ggplot2::element_rect(fill = "grey15"),
      panel.grid.major  = ggplot2::element_line(color = "grey25")
    )
}
