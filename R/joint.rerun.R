#' Rerun a model fitted with INLAjoint
#'
#' @description Reruns inla() for an object of class \code{INLAjoint} returned by the \code{joint}
#'   function. The rerun starts with posterior values from previous run and can sometimes improve the model fit
#'   (for very complex models or unstable parameter estimates due to low information in the data)
#' @param model an object containing a model fitted with the joint() function.
#' @param ... Extra arguments.
#'
#' @seealso \code{\link{joint}}.
#' @return An object of class \code{INLAjoint} containing a model fitted with the joint() function.
#'
#' @export

joint.rerun <- function(model, ...){
  rerun.model <- INLA::inla.rerun(model)
  CLEANoutput <- c('summary.lincomb','mfarginals.lincomb','size.lincomb',
                   'summary.lincomb.derived','marginals.lincomb.derived','size.lincomb.derived','offset.linear.predictor',
                   'model.spde2.blc','summary.spde2.blc','marginals.spde2.blc','size.spde2.blc','model.spde3.blc','summary.spde3.blc',
                   'marginals.spde3.blc','size.spde3.blc','logfile','Q','graph','ok','model.matrix')
  rerun.model[CLEANoutput] <- NULL
  rerun.model$cureVar <- model$cureVar
  rerun.model$variant <- model$variant
  rerun.model$famLongi <- model$famLongi
  rerun.model$corLong <- model$corLong
  rerun.model$longOutcome <- model$longOutcome
  rerun.model$survOutcome <- model$survOutcome
  rerun.model$id <- model$id
  rerun.model$timeVar <- model$timeVar
  rerun.model$REstruc <- model$REstruc
  rerun.model$REstrucS <- model$REstrucS
  rerun.model$basRisk <- model$basRisk
  rerun.model$call <- model$call
  rerun.model$dataLong <- model$dataLong
  rerun.model$dataSurv <- model$dataSurv
  class(rerun.model) <- c("INLAjoint", "inla")
  return(rerun.model)
}
