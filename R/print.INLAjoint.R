#' @export

print.INLAjoint <- function(x, sdcor=FALSE, hazr=FALSE, ...){
  summary(x,sdcor=sdcor, hazr=hazr)
}
