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
#' @param ch.terms a vector of column names that should be used as conditional heteroscedasticity covariates
#' @return A vector containing the outputated regression coefficient estimates
#' @export
mixedtobit <- function(formula, data, M, left = -1, id, ch.terms = NULL)
{
  # Note: right now, only supports
  # (a) full mixed effects model (so, fixed AND random effect for each parameter provided)
  # (b) left-censoring only
  # TODO: Remove these limitations, especially if we are to submit this to CRAN.
  #       These aren't terribly difficult, just time-consuming to implement.

  X <- model.matrix(formula, data = data)
  y.name <- all.vars(formula)[1]

  fn <- function(fn.data)
  {
    fn.X <- model.matrix(formula, data = fn.data)

    if(!is.null(ch.terms))
    {
      fn.formula <- fn.data[, y.name] ~ -1 + fn.X | as.matrix(fn.data[, ch.terms])
    } else
    {
      fn.formula <- fn.data[, y.name] ~ -1 + fn.X | 1
    }

    m <- crch(fn.formula,
              link.scale = "quadratic", left = left,
              data = fn.data)

    beta <- m$coefficients$location
    names(beta) <- colnames(X)

    Sigma <- m$vcov[1:ncol(fn.X), 1:ncol(fn.X)]

    ch.var <- m$coefficients$scale
    names(ch.var)[-1] <- ch.terms

    return(list(beta = beta, Sigma = Sigma,
                ch.var = ch.var))
  }


  mo <- try(multiout(fn = fn, M = M, data = data, id = id, leave.as.list = TRUE))
  tries <- 1
  while(class(mo) == "try-error" && tries < 10)
  {
    mo <- try(multiout(fn = fn, M = M, data = data, id = id, leave.as.list = TRUE))
  }

  if(tries >= 10)
  {
    stop("Too many failed attempts")
  }

  betas <- Reduce(function(a,b){rbind(a, b$beta)}, x = mo[-1], init = mo[[1]]$beta)
  beta.hat <- colMeans(betas)

  Sigma.sum <- Reduce(function(a, b){a + b$Sigma}, x = mo[-1], init = mo[[1]]$Sigma)

  beta.var <- var(betas)

  Sigma.of.est <- Sigma.sum / M - beta.var # See (Follman, 2003)

  ch.var.sum <- Reduce(function(a, b){a + b$ch.var}, x = mo[-1], init = mo[[1]]$ch.var)
  ch.var.hat <- ch.var.sum / M

  return(list(beta = beta.hat,
              Sigma = Sigma.of.est,
              ch.var = ch.var.hat,
              mean.Sigmas = Sigma.sum / M,
              beta.var = beta.var))
}
