#' @export



print.summary.INLAjoint <- function(x, ...){
  if (!class(x)=="summary.INLAjoint"){
    stop("Please provide an object of class 'summary.INLAjoint' (obtained with joint() function).\n")
  }
  Nerr <- 1 #  identify error terms
  for(i in 1:x$NLongi){
    if(x$NLongi==1){
      cat(paste0("Longitudinal outcome (", x$famLongi[i], ")\n"))
    }else{
      if(i==1){
        cat(paste0("Longitudinal outcome (L", i, ", ", x$famLongi[i],")\n"))
      }else{
        cat(paste0("\nLongitudinal outcome (L", i, ", ", x$famLongi[i],")\n"))
      }
    }
    print(round(x$FixedEff[[i]], 4))
    if(x$NLongi==1){
      if(!x$sdcor) cat("\nRandom effect variance \n") else cat("\nRandom effect standard deviation\n")
      print(round(x$ReffList[[i]], 4))
    }else if(x$NRand==x$NLongi){
      if(!x$sdcor) cat(paste0("\nRandom effects variance-covariance (L", i, ")\n")) else cat(paste0("\nRandom effects standard deviation / correlation (L", i, ")\n"))
      print(round(x$ReffList[[i]], 4))
    }else if(i==x$NLongi){
      if(!x$sdcor) cat(paste0("\nRandom effects variance-covariance\n")) else cat(paste0("\nRandom effects standard deviation / correlation\n"))
      print(round(x$ReffList[[1]], 4))
    }
  }
  if(x$NSurv>0){
    for(i in 1:x$NSurv){
      if(x$NSurv==1){
        cat("\nSurvival outcome\n")
      }else{
        cat(paste0("\nSurvival outcome (S", i, ")\n"))
      }
      print(round(x$SurvEff[[i]], 4))
    }
  }
  if(dim(x$AssocLS)[1]>0){
    cat("\nAssociation longitudinal - survival\n")
    print(round(x$AssocLS, 4))
  }
  cat("\n")
  print(x$mlik[,1])
  if(!is.null(x$dic)) cat(paste0("\nDeviance Information Criterion: "), x$dic)
  if(!is.null(x$waic)) cat(paste0("\nWidely applicable Bayesian information criterion: "), x$waic)
  cat(paste0("\nComputation time: ", round(x$cpu.used[4], 2), " seconds"))
  cat("\n")
  invisible(x)
}