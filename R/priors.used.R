#' Get prior distributions used for an object of class INLAjoint
#'
#' @description Returns the prior distribution for all the parameters of an object of class \code{INLAjoint} returned by the \code{joint}
#'   function that fits a joint model to multivariate longitudinal and
#'   time-to-event data.
#'
#' @seealso \code{\link{joint}}, \code{\link{summary.INLAjoint}}.
#' @return A list of the prior distribution used for a model of class INLAjoint fitted with joint().
#'
#' @export

priors.used <- function(x){
  if(!("INLAjoint" %in% class(x))) stop("The function priors.used only applies to an object of class 'INLAjoint' as returned by the 'joint' function of the 'INLAjoint' package.")
  class(x) <- "inla"
  print(inla.priors.used(x))
}


