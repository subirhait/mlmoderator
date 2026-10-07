# mlmoderator (development version)

* `plot.mlm_variance_decomp()`: the subtitle labelled the random-slope
  standard deviation as tau11, which is the random-slope variance. It now
  states both quantities explicitly ("random-slope SD = ..., variance = ...").
  The subtitle is also ASCII only: the Unicode tau it used was dropped by the
  default `pdf()` graphics device and printed as "..". Plotted intervals and
  all numerical results are unchanged.
* Tests: the shared test fixture no longer writes its data to the global
  environment. The data are embedded in the model call and read back with
  `lme4::getData()`, so the tests run under current testthat, which isolates
  tests from the global environment.

# mlmoderator 0.3.0

## Corrections to inference (results change)

* **Degrees of freedom.** All tests and intervals now use Satterthwaite
  degrees of freedom by default (via 'lmerTest'). Versions up to 0.2.1 used
  N - p (level-1 rows minus fixed effects), which treats observations as
  independent and is anti-conservative for cross-level interactions: in a
  simulation with 10 clusters and no true interaction it rejected at about
  twice the nominal rate. The new `df_method` argument offers
  `"satterthwaite"` (default), `"kenward-roger"`, `"between"`
  (J - q - 1; Snijders and Bosker, 2012), and `"residual"` (the old rule,
  kept only for reproducing earlier results).
* **Johnson-Neyman boundaries** are now exact: closed-form roots of the
  Bauer-Curran quadratic for constant-df methods, and root-finding on the
  exact criterion for Satterthwaite and Kenward-Roger. Previously they were
  interpolated from p-values on a grid. `mlm_jn()` also returns
  `jn_bounds_all` (all real roots, including those outside the data) for
  constant-df methods.
* **Confidence bands in `mlm_plot()`** now use the full fixed-effects design
  row. Previously covariates other than the predictor and moderator were
  omitted from the standard error of the fitted values.
* **Moderator probe values** for a cluster-level moderator are now the mean
  and SD across clusters, not across observations (new `modx.level`
  argument; `"observation"` restores the old behaviour).

## `mlm_sensitivity()` (breaking)

* The ICC-shift analysis and `robustness_index` were removed. The
  design-effect rescaling it used applies to means under a random-intercept
  model, not to cross-level interactions, and the adjusted Johnson-Neyman
  boundary it reported was not correct. `icc_range` and `icc_grid` now give
  a warning and are ignored; `loco = FALSE` is an error.
* Leave-one-cluster-out refits now use `update()` on the original data, so
  they keep the original REML/ML setting, weights, and formula. Previously
  they always used ML (so every "change" included the ML/REML difference)
  and failed silently for formulas with transformed variables.
* The influence measure is now DFBETA, the change in the interaction divided
  by its full-data SE, flagged at 2/sqrt(J). The previous
  "Cook's-distance-style" measure divided by the SD of the changes and so
  flagged clusters by construction. Failed refits are reported with their
  error messages.

## `mlm_variance_decomp()` (breaking)

* `pct_random` was removed: it divided a population variance by a sampling
  variance and so grew with sample size. The output now reports `tau11`
  (variance), `tau11_sd`, and `df` explicitly, and the documentation
  describes the result as confidence intervals for the average slope versus
  prediction intervals for a new cluster.
* The random slope is taken from the grouping factor that carries `pred`.

## Other

* Tests now cover every exported function, including checks against
  'lmerTest', 'emmeans', and 'pbkrtest'.
* New vignette with the public High School and Beyond data ('mlmRev').
* Removed a reference to a non-existent vignette.
* URLs point to github.com/subirhait/mlmoderator.

# mlmoderator 0.2.1

* Previous CRAN release.
