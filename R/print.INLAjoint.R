#' @export

print.INLAjoint <- function (x, sdcor=FALSE, ...) {
  summary(x,sdcor=sdcor)
  invisible(x)
}