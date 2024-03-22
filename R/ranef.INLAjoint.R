#' Extracts random effects values from a given model fitted with INLAjoint
#'
#' @description This function extracts random effects values from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export
#' @importFrom nlme ranef

ranef.INLAjoint <- function(object, ...){
  arguments <- list(...)
  if(is.null(arguments$hr)) hr=F else hr=arguments$hr
  if(!is.null(arguments$hazr)) hr=arguments$hazr
  if(is.null(arguments$sdcor)) sdcor=F else sdcor=arguments$sdcor
  SUM <- summary(object, hr=hr, sdcor=sdcor)
  res <- NULL
  if(inherits(SUM$ReffList, "list")){
    if(length(SUM$ReffList)>0){
      for(i in 1:length(SUM$ReffList)){
        res <- rbind(res, SUM$ReffList[[i]])
      }
    }
  }else if(inherits(SUM$ReffList, "data.frame")){
    if(dim(SUM$ReffList)[1]>0){
      res <- rbind(res, SUM$ReffList)
    }
  }
  return(res)
}








