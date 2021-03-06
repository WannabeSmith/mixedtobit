% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixedtobit.R
\name{mixedtobit}
\alias{mixedtobit}
\title{Estimate tobit regression model parameters with mixed effects}
\usage{
mixedtobit(formula, data, M, left = -1, id, ch.terms = NULL,
  beta.only = FALSE)
}
\arguments{
\item{formula}{a regression formula describing the relationship between the response and the predictors}

\item{data}{a data.frame containing the response and predictors of interest, as well as cluster IDs}

\item{M}{the (positive integer) number of outputations (resamples) to perform}

\item{left}{the location of left censoring}

\item{id}{a string of the column name in "data" specifying the cluster ID}

\item{ch.terms}{a vector of column names that should be used as conditional heteroscedasticity covariates}

\item{beta.only}{a logical variable indicating whether or not to include just beta (this is useful for simulations with memory requirements)}
}
\value{
A vector containing the outputated regression coefficient estimates
}
\description{
Use multiple outputation (a.k.a. within-cluster resampling) to estimate tobit regression model parameters with mixed effects
}
