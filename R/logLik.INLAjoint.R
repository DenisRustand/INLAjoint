#' Extracts log-likelihood value from a given model fitted with INLAjoint
#'
#' @description This function extracts log-likelihood value from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

logLik.INLAjoint <- function(object, ...){
  arguments <- list(...)
  return(object$mlik)
}








