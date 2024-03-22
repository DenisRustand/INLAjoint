#' Extracts fixed effects values from a given model fitted with INLAjoint
#'
#' @description This function extracts fixed effects values from INLAjoint objects.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param ... Extra arguments.
#'
#' @export
#' @importFrom nlme fixef

fixef.INLAjoint <- function(object, ...){
  arguments <- list(...)
  if(is.null(arguments$hr)) hr=F else hr=arguments$hr
  if(!is.null(arguments$hazr)) hr=arguments$hazr
  if(is.null(arguments$sdcor)) sdcor=F else sdcor=arguments$sdcor
  SUM <- summary(object, hr=hr, sdcor=sdcor)
  res <- NULL
  if(inherits(SUM$FixedEff, "list")){
    if(length(SUM$FixedEff)>0){
      for(i in 1:length(SUM$FixedEff)){
        if(length(grep("Res. err.", rownames(SUM$FixedEff[[i]])))>0){
          res <- rbind(res, SUM$FixedEff[[i]][-grep("Res. err.", rownames(SUM$FixedEff[[i]])),])
        }else{
          res <- rbind(res, SUM$FixedEff[[i]])
        }
      }
    }
  }else if(inherits(SUM$FixedEff, "data.frame")){
    if(dim(SUM$FixedEff)[1]>0){
      res <- rbind(res, SUM$FixedEff)
    }
  }
  if(inherits(SUM$SurvEff, "list")){
    if(length(SUM$SurvEff)>0){
      for(i in 1:length(SUM$SurvEff)){
        if(length(grep("Baseline", rownames(SUM$SurvEff[[i]])))>0){
          res <- rbind(res, SUM$SurvEff[[i]][-grep("Baseline", rownames(SUM$SurvEff[[i]])),])
        }else if(length(grep("Exponential", rownames(SUM$SurvEff[[i]])))>0){
          res <- rbind(res, SUM$SurvEff[[i]][-grep("Exponential", rownames(SUM$SurvEff[[i]])),])
        }else if(length(grep("Weibull", rownames(SUM$SurvEff[[i]])))>0){
          res <- rbind(res, SUM$SurvEff[[i]][-grep("Weibull", rownames(SUM$SurvEff[[i]])),])
        }else{
          res <- rbind(res, SUM$SurvEff[[i]])
        }
      }
    }
  }else if(inherits(SUM$SurvEff, "data.frame")){
    if(dim(SUM$SurvEff)[1]>0){
      res <- rbind(res, SUM$SurvEff)
    }
  }
  if(dim(SUM$AssocSS)[1]) res <- rbind(res, SUM$AssocSS)
  if(dim(SUM$AssocLS)[1]) res <- rbind(res, SUM$AssocLS)
  return(res)
}








