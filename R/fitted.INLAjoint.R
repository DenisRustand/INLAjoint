#' Extracts fitted values from a given model fitted with INLAjoint
#'
#' @description This function extracts fitted values from INLAjoint objects.
#' Values are associated to a name to keep track of the outcome related to each fitted value.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export

fitted.INLAjoint <- function(object, ...){
  arguments <- list(...)
  PRED <- object$summary.fitted.values$mean
    OUtc <- as.data.frame(object$.args$data$Yjoint)
    names(PRED) <- sapply(1:length(PRED), function(x) colnames(OUtc)[which(!is.na(OUtc[x,]))])
    return(PRED)
}








