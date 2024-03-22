#' Extracts model coefficients from a given model fitted with INLAjoint
#'
#' @description This function extracts model coefficients from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

coef.INLAjoint <- function(object, ...){
  arguments <- list(...)
  FEF <- object$summary.fixed
  HYP <- object$summary.hyperpar
  OUtc <- list("Fixed_effects" = FEF, "Hyperparameters" = HYP)
  return(OUtc)
}








