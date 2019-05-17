#' Estimate tobit regression model parameters with mixed effects
#'
#' Use multiple outputation (a.k.a. within-cluster resampling) to estimate tobit regression model parameters with mixed effects
#'
#' @importFrom multiout multiout
#' @importFrom crch crch
#' @importFrom stats model.matrix
#' @param formula a regression formula describing the relationship between the response and the predictors
#' @param data a data.frame containing the response and predictors of interest, as well as cluster IDs
#' @param M the (positive integer) number of outputations (resamples) to perform
#' @param left the location of left censoring
#' @param id a string of the column name in "data" specifying the cluster ID
#' @return A vector containing the outputated regression coefficient estimates
#' @export
mixedtobit <- function(formula, data, M, left = -1, id)
{
  # Note: right now, only supports
  # (a) full mixed effects model (so, fixed AND random effect for each parameter provided)
  # (b) same variance for each random effect.
  # (c) left-censoring only
  # TODO: Remove these limitations, especially if we are to submit this to CRAN.
  #       These aren't terribly difficult, just time-consuming to implement.

  X <- model.matrix(formula, data = data)
  data$xTx <- rowSums(X)

  y.name <- all.vars(formula)[1]

  fn <- function(fn.data)
  {
    fn.X <- model.matrix(formula, data = fn.data)

    m <- crch(fn.data[, y.name] ~ -1 + fn.X | xTx,
              link.scale = "quadratic", left = -1,
              data = fn.data)

    beta <- m$coefficients$location
    names(beta) <- colnames(X)

    return(beta)
  }

  beta.ests <- multiout(fn = fn, M = M, data = data, id = id)

  return(beta.ests)
}
