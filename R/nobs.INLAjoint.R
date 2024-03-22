#' Extracts number of observations of each composant from a given model fitted with INLAjoint
#'
#' @description This function extracts number of observations of each composant from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

nobs.INLAjoint <- function(object, ...){
  arguments <- list(...)
  OUtc <- as.data.frame(object$.args$data$Yjoint)
  return(table(sapply(1:length(object$summary.fitted.values$mean), function(x) colnames(OUtc)[which(!is.na(OUtc[x,]))])))
}








