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
#' print by specifying the names or corresponding letters in the argument 'which'.
#'
#' @param x an object with the output of the the \link{joint} function
#' @param which name of required plots. Default is "all".
#'  It can be a character or named list.
#'  If it is a list, each element can be either character or numeric to
#'  select from the elements to be visualized.
#'  Ex.: which = list(Outcomes = "L1") and list(Outcomes = 1)
#'  will produce the same output.
#' @param ... Extra arguments.
#' @import ggplot2
#' @importFrom grDevices dev.new
#' @export
print.plot.INLAjoint <- function(x,
                                 which = c("all", "Outcomes", "Covariances",
                                           "Associations", "Baseline"), ...) {
  arguments <- list(...)
  stopifnot(class(which) %in% c("character", "list"))
  if(dev.interactive()){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  wch <- c("all", "Outcomes", "Covariances", "Associations", "Baseline")
  if(is.character(which)) {
    which <- unique(which)
    for(i in 1:length(which)) {
      which[i] <- match.arg(which[i], wch)
    }
    if(any(which=="all")) {
      iwhich <- vector("list", length(x))
      iw1 <- 1:length(x)
      for(i in iw1) {
        if(all(class(x[[i]]) == "list")) {
          iwhich[[i]] <- 1:length(x[[i]])
        }
      }
    } else {
      iw1 <- pmatch(which, wch[-1])
      iw1 <- iw1[!is.na(iw1)]
      if(length(iw1) == 0) {
        stop('Invalid "which" argument!')
      }
      iwhich <- vector("list", length(iw1))
      for(i in 1:length(iw1)) {
        if(all(class(x[[i]]) == "list")) {
          iwhich[[i]] <- 1:length(x[[iw1[i]]])
        }
      }
    }
  }
  if(is.list(which)) {
    iw1 <- pmatch(names(which), wch[-1])
    if(all(is.na(iw1))) {
      stop('Invalid "names(which)"')
    }
    if(any(is.na(iw1))) {
      warning('Some elements of "names(which)" are invalid"!')
    }
    iw1 <- iw1[!is.na(iw1)]
    iwhich <- vector("list", length(iw1))
    for (i in 1:length(which)) {
      if(all(class(x[[i]]) == "list")) {
        if(is.character(which[[i]])) {
          ij <- pmatch(which[[i]], names(x[[iw1[i]]]))
          if(all(is.na(ij))) {
            warning('Found invalid element(s) in the element', i, 'of "which"!')
          }
          iwhich[[i]] <- ij[!is.na(ij)]
        } else {
          if(!is.numeric(which[[i]])) {
              warning('Element', i, 'of "which" is invalid!')
          } else {
            iwhich[[i]] <- intersect(which[[i]], 1:length(x[[i]]))
            if(length(iwhich[[i]])<length(which[[i]])) {
                warning('Found invalid element(s) in the element', i, 'of "which"!')
            }
            if(length(iwhich[[i]])==0) {
               warning('Element', i, 'of "which" is invalid!')
            }
          }
        }
      }
    }
  }
  for(i in iw1) {
    if(all(class(x[[i]]) == "list")) {
      for(j in iwhich[[i]]) {
        print(x[[i]][[j]])
        dev.flush()
      }
    } else {
      print(x[[i]])
    }
  }
  devAskNewPage(FALSE)
}
