#' INLAjoint
#'
#' @description INLAjoint is a package that fits joint models for multivariate longitudinal markers
#' (with various distributions available) and survival outcomes (possibly accounting
#' for competing risks) with Integrated Nested Laplace Approximations (INLA). The
#' flexible and user friendly function joint() facilitates the use of the fast and
#' reliable inference technique implemented in INLA package for joint modeling. More
#' details are given in the help page of the joint function (accessible via ?joint in
#' the R console), the vignette associated to the joint() function (accessible via
#' vignette("INLAjoint") in the R console).
#'
#' Contact: \email{INLAjoint@gmail.com}
#'
#' @docType package
#' @name INLAjoint
INLAjoint <- function() {
  message("Welcome to the INLAjoint package!")
  #utils::browseVignettes("INLAjoint")
}
