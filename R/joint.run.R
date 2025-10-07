#' Run a model fitted with INLAjoint
#'
#' @description Runs inla() for an object of class \code{INLAjoint} returned by the \code{joint}
#'   function with argument `run` set to FALSE. The rerun starts with posterior values from
#'   previous run and can sometimes improve the model fit
#'   (for very complex models or unstable parameter estimates due to low information in the data)
#' @param model an object containing a model fitted with the joint() function.
#' @param silentMode boolean to display messages about the fit procedure. Default is FALSE.
#' @param class defines the class of the object created. Default is "INLAjoint" but can be switched to "inla".
#' @param ... Extra arguments.
#'
#' @seealso \code{\link{joint}}.
#' @return An object of class \code{INLAjoint} containing a model fitted with the joint() function.
#'
#' @export

joint.run <- function(model, silentMode=FALSE, class="INLAjoint", ...){
  if(model$run) stop("Model already did run!")
  class(model) <- "inla"
  if(!silentMode) message("Fit model...")
  res <- INLA::inla.rerun(model, plain=TRUE)
  if(model$control$rerun){
    CT1 <- res$cpu.used[4]
    res <- INLA::inla.rerun(res)
    res$cpu.used[4] <- res$cpu.used[4] + CT1 # account for first fit in total computation time
  }
 # not compatible with NL effects yet (could be aded upon request, needs to do the stepwise procedure here then...)
  if(!is.null(model$assoc_Names)){
    for(a_s in model$assoc_Names){
      if(length(grep("SRE_ind", a_s))==0){
        res$dic$local.dic[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- 0
        res$dic$local.p.eff[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- 0
        res$waic$local.waic[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- 0
        res$waic$local.p.eff[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- 0
        res$dic$dic <- sum(na.omit(res$dic$local.dic))
        res$dic$p.eff <- sum(na.omit(res$dic$local.p.eff))
        res$waic$waic <- sum(na.omit(res$waic$local.waic))
        res$waic$p.eff <- sum(na.omit(res$waic$local.p.eff))
        if(length(res$cpo$cpo)>0){
          res$cpo$cpo[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- NA
        }
        if(length(res$cpo$pit)>0){
          res$cpo$pit[which(!is.na(res$.args$data$Yjoint[[a_s]]))] <- NA
        }
      }
    }
  }
  if(length(res$misc$warnings)>0 & "Skewne" %in% substr(res$misc$warnings, 1, 6)) warning("The hyperparameters skewness correction seems abnormal, this can be a sign of an ill-defined model and/or issues with the fit.")
  if(length(res$misc$warnings)>0 & "Stupid" %in% substr(res$misc$warnings, 1, 6)) warning("Stupid local search strategy used: This can be a sign of a ill-defined model and/or non-informative data.")
  if(TRUE %in% c(abs(res$misc$cor.intern[upper.tri(res$misc$cor.intern)])>0.99))
    warning("Internal correlation between hyperparameters is abnormally high, this is a sign of identifiability issues / ill-defined model. ")
  CLEANoutput <- c('offset.linear.predictor',
                   # 'summary.lincomb','mfarginals.lincomb','size.lincomb', 'summary.lincomb.derived','marginals.lincomb.derived','size.lincomb.derived',
                   'model.spde2.blc','summary.spde2.blc','marginals.spde2.blc','size.spde2.blc','model.spde3.blc','summary.spde3.blc',
                   'marginals.spde3.blc','size.spde3.blc','Q','graph','ok','model.matrix')
  res[CLEANoutput] <- NULL
  res$run <- TRUE # the model did run
  # if(!is.null(res$lincomb)) res$lincomb <- NULL
  if(model$is_Surv) res$cureVar <- model$cureVar
  if(model$is_Surv) res$variant <- model$variant
  if(model$is_Surv) res$cutpoints <- model$cutpoints
  if(model$is_Surv) res$NbasRisk <- model$NbasRisk
  if(model$is_Surv & !is.null(model$SurvInfo)) res$SurvInfo <- model$SurvInfo
  if(model$is_Long) res$famLongi <- model$famLongi
  if(model$is_Long) res$corLong <- model$corLong else res$corLong <- FALSE
  if(model$is_Long) res$control.link <- model$control.link
  if(model$is_Long) res$longOutcome <- model$longOutcome
  if(model$is_Surv) res$survOutcome <- model$survOutcome
  if(!is.null(model$formulaAssocInfo)) res$assoc <- model$assoc
  if(!is.null(model$assoc_Names)) res$assoc_Names <- model$assoc_Names
  res$id <- model$id
  res$timeVar <- model$timeVar
  if(!is.null(model$range)) res$range <- model$range
  if(!is.null(model[["REstruc"]])) res$REstruc <- model[["REstruc"]]
  if(!is.null(model$REstruc)) res$mat_k <- model$mat_k
  if(!is.null(model$control$fixRE)) res$fixRE <- model$fixRE
  if(!is.null(model$control$strata)) res$strata <- model$strata
  if(!is.null(model$lonFacChar)) res$lonFacChar <- model$lonFacChar
  if(!is.null(model$survFacChar)) res$survFacChar <- model$survFacChar
  if(!is.null(model[["REstruc"]])) res$corRE <- model$corRE # switch for diagonal/correlated random effects within a longitudinal marker
  if(!is.null(model$REstrucS)) res$REstrucS <- model$REstrucS
  res$formSurv <- model$formSurv
  res$formLong <- model$formLong
  if(!is.null(model$NLcov_name)) res$NLinfo <- model$NLinfo
  res$basRisk <- model$basRisk
  res$priors_used <- model$priors_used
  res$call <- model$call
  res$dataLong <- model$dataLong
  res$dataSurv <- model$dataSurv
  if(!silentMode) message("...done!")
  if(class=="INLAjoint") class(res) <- c("INLAjoint", "inla") else class(res) <- "inla"
  return(res)
}
