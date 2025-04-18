% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint.run.R
\name{joint.run}
\alias{joint.run}
\title{Run a model fitted with INLAjoint}
\usage{
joint.run(model, silentMode = FALSE, class = "INLAjoint", ...)
}
\arguments{
\item{model}{an object containing a model fitted with the joint() function.}

\item{silentMode}{boolean to display messages about the fit procedure. Default is FALSE.}

\item{class}{defines the class of the object created. Default is "INLAjoint" but can be switched to "inla".}

\item{...}{Extra arguments.}
}
\value{
An object of class \code{INLAjoint} containing a model fitted with the joint() function.
}
\description{
Runs inla() for an object of class \code{INLAjoint} returned by the \code{joint}
  function with argument `run` set to FALSE. The rerun starts with posterior values from
  previous run and can sometimes improve the model fit
  (for very complex models or unstable parameter estimates due to low information in the data)
}
\seealso{
\code{\link{joint}}.
}
