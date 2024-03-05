#' Extracts formula from a given model fitted with INLAjoint
#'
#' @description This function extracts formula from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.

#' @export

formula.INLAjoint <- function(object, ...){
  OUtc <- list("Formula.INLAjoint" = list("Longitudinal" = object$formLong, "Survival" = object$formSurv), "Formula.INLA" = object$.args$formula)
  return(OUtc)
}








