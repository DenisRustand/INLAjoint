#' Prints plot the output from a multivariate joint model for longitudinal and/or survival data
#'
#' @description The plots of a model are grouped by categories, first are the fixed effects and
#' residual error of longitudinal and survival submodels, referred to as 'Outcomes' or 'O'. Then the
#' variance-covariance of random effects (or standard deviations and correlations when the
#' argument 'sdcor' is set to TRUE in the call of the plot function), referred to as
#' 'Covariances' or 'C'. Association parameters referred to as 'Associations' or 'A' for
#' linear associations and 'NL_Associations' or 'N' for non-linear associations. Baseline
#' hazard curves referred to as 'Baseline' or 'B' and baseline hazard related parameters
#' referred to as 'ParamBaseline' or 'P'. It is possible to select specific plots to
#' print by specifying the names or corresponding lettes in the argument 'which'.
#'
#' @param x an object with the output of the the \link{joint} function
#' @param which name of required plots. Default is all.
#'
#' @import ggplot2
#' @importFrom grDevices dev.new
#' @export
print.plot.INLAjoint <- function(x, which="all", ...) {
  if(which=="all"){
    if(dev.interactive()){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    for(i in 1:length(x)) {
      # if(i>1) dev.new()
      if(all(class(x[[i]]) == "list")) {
        for(j in 1:length(x[[i]])) {
          print(x[[i]][[j]])
          dev.flush()
        }
      } else {
        print(x[[i]])
      }
    }
  }else{
    print("WIP")
  }
  devAskNewPage(FALSE)
}
