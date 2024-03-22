#' Extracts family from a given model fitted with INLAjoint
#'
#' @description This function extracts family from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

family.INLAjoint <- function(object, ...){
  arguments <- list(...)
  if(!is.null(object$formLong)){
    Long <- object$famLongi
  }else{
    Long <- NULL
  }
  if(!is.null(object$formSurv)){
    Surv <- unlist(object$basRisk)
  }else{
    Surv <- NULL
  }
  return(list("Longitudinal" = Long, "Survival" = Surv))
}








