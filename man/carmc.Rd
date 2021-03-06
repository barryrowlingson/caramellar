% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carmc.R
\name{carmc}
\alias{carmc}
\title{Conditional Autoregressive Space-time MCMC}
\usage{
carmc(cases, space_time, linear_model, neighbours, mc_control,
  mcmc = mcmc.defaults(), note = NA, progress)
}
\arguments{
\item{cases}{data frame of case data}

\item{space_time}{formula of location_id~date}

\item{linear_model}{formula of caseflag~covariates}

\item{neighbours}{neighbourhood structure}

\item{mc_control}{control mcmc loops, thinning, burn-in}

\item{mcmc}{priors, proposals, etc}

\item{note}{a text note to add to the document}

\item{progress}{passed to txtProgressBar to create a progress bar}
}
\value{
an object of class "carmc"
}
\description{
Space-time Conditional Autoregressive MCMC
}
\details{
Fit the space-time conditional AR model via MCMC
}
\author{
Emanuele Giorgi, Barry Rowlingson, Alison Hale
}

