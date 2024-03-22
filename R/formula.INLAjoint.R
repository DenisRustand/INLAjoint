#' Extracts formula from a given model fitted with INLAjoint
#'
#' @description This function extracts formula from INLAjoint objects.
#'
#' @param x an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

formula.INLAjoint <- function(x, ...){
  arguments <- list(...)
  OUtc <- list("Formula.INLAjoint" = list("Longitudinal" = x$formLong, "Survival" = x$formSurv), "Formula.INLA" = x$.args$formula)
  return(OUtc)
}








