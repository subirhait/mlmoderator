# mlmoderator

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mlmoderator)](https://CRAN.R-project.org/package=mlmoderator)
<!-- badges: end -->

`mlmoderator` probes, plots, and checks cross-level interactions in
two-level models fitted with `lme4::lmer()`.

| Function | What it does |
|---|---|
| `mlm_center()` | Grand-mean, group-mean, or within/between centring |
| `mlm_probe()` | Simple slopes at chosen moderator values |
| `mlm_jn()` | Johnson-Neyman boundaries (closed form or exact root-finding) |
| `mlm_plot()` | Interaction plot with confidence bands |
| `mlm_surface()` | Contour plot of predicted outcomes over predictor x moderator |
| `mlm_summary()` | Interaction test, simple slopes, and JN region together |
| `mlm_variance_decomp()` | Confidence intervals for the average slope vs. prediction intervals for a new cluster |
| `mlm_sensitivity()` | Leave-one-cluster-out influence (DFBETA) on the interaction |

## Degrees of freedom

Inference for a cross-level interaction draws its information from the
clusters. All tests and intervals therefore use Satterthwaite degrees of
freedom by default (via **lmerTest**), with Kenward-Roger (`"kenward-roger"`)
and a between-cluster rule (`"between"`, J - q - 1) as alternatives, through
the `df_method` argument. Versions before 0.3.0 used N - p, which is
anti-conservative with few clusters; it remains available as
`df_method = "residual"` for reproducing earlier results.

## Installation

``` r
install.packages("mlmoderator")
```

## Example

``` r
library(mlmoderator)
library(lme4)

data(school_data)
mod <- lmer(math ~ ses * climate + gender + (1 + ses | school),
            data = school_data)

mlm_summary(mod, pred = "ses", modx = "climate")
mlm_plot(mod, pred = "ses", modx = "climate")
plot(mlm_jn(mod, pred = "ses", modx = "climate"))
mlm_variance_decomp(mod, pred = "ses", modx = "climate")
mlm_sensitivity(mod, pred = "ses", modx = "climate")
```

`vignette("hsb-workflow")` works through a full analysis of the public High
School and Beyond data.

## Scope

The package describes and checks the fitted model. It does not address
unmeasured confounding of the interaction, and it currently supports
Gaussian `lmer()` models with the cluster defined by the first grouping
factor.
