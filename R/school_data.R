#' Simulated school achievement dataset
#'
#' A simulated two-level dataset with students nested within schools, designed
#' to illustrate multilevel moderation analysis. The true data-generating model
#' includes a cross-level interaction between student socioeconomic status and
#' school climate.
#'
#' @format A data frame with 3,000 rows and 6 variables:
#' \describe{
#'   \item{school}{A factor indicating the school identifier (1--100).}
#'   \item{student}{An integer indicating the student identifier (1--3000).}
#'   \item{math}{A numeric mathematics achievement score.}
#'   \item{ses}{A numeric student socioeconomic status variable (standardized).}
#'   \item{climate}{A numeric school climate rating (standardized level-2 variable).}
#'   \item{gender}{A factor indicating student gender with levels `"female"` and `"male"`.}
#' }
#'
#' @details
#' The data were generated from the model
#' \deqn{math_{ij} = 50 + 1.5\,ses_{ij} + 0.8\,climate_j +
#' 0.5\,ses_{ij} climate_j + u_{0j} + u_{1j} ses_{ij} + e_{ij}.}
#'
#' The level-2 random effects were generated as \eqn{u_{0j} \sim N(0, 9)} and
#' \eqn{u_{1j} \sim N(0, 0.25)}, and the level-1 residuals were generated as
#' \eqn{e_{ij} \sim N(0, 25)}.
#'
#' @examples
#' data(school_data)
#' head(school_data)
#' str(school_data)
#'
#' \donttest{
#' library(lme4)
#' mod <- lmer(math ~ ses * climate + gender + (1 + ses | school),
#'             data = school_data)
#' mlm_probe(mod, pred = "ses", modx = "climate")
#' }
"school_data"
