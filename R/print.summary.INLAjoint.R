#' @export



print.summary.INLAjoint <- function(x, ...){
  if (!inherits(x, "summary.INLAjoint")){
    stop("Please provide an object of class 'summary.INLAjoint' (obtained with joint() function).\n")
  }
  if(x$NLongi>0){
    Nerr <- 1 #  identify error terms
    for(i in 1:x$NLongi){
      if(x$NLongi==1){
        if(x$NSurv==0) rownames(x$FixedEff[[1]]) <- sapply(rownames(x$FixedEff[[1]]), function(x) gsub("_L1", "", x))
        cat(paste0("Longitudinal outcome (", x$famLongi[i], ")\n"))
      }else{
        if(i==1){
          cat(paste0("Longitudinal outcome (L", i, ", ", x$famLongi[i],")\n"))
        }else{
          cat(paste0("\nLongitudinal outcome (L", i, ", ", x$famLongi[i],")\n"))
        }
      }
      print(round(x$FixedEff[[i]], 4))
      if(x$NLongi==1 & x$NSurv==0) rownames(x$ReffList[[1]]) <- sapply(rownames(x$ReffList[[1]]), function(x) gsub("_L1", "", x))
      if(x$NRand==x$NLongi){
        if(!x$sdcor) cat(paste0("\nRandom effects variance-covariance (L", i, ")\n")) else cat(paste0("\nRandom effects standard deviation / correlation (L", i, ")\n"))
        print(round(x$ReffList[[i]], 4))
      }else if(i==x$NLongi){
        if(!x$sdcor) cat(paste0("\nRandom effects variance-covariance\n")) else cat(paste0("\nRandom effects standard deviation / correlation\n"))
        print(round(x$ReffList[[1]], 4))
      }
    }
  }else if(!is.null(x[["ReffList"]])){ # in case of random effects only in longitudinal parts
    if(!x$sdcor) cat(paste0("\nRandom effects variance-covariance\n")) else cat(paste0("\nRandom effects standard deviation / correlation\n"))
    print(round(x$ReffList[[1]], 4))
  }
  if(x$NSurv>0){
    for(i in 1:x$NSurv){
      if(x$NSurv==1){
        if(x$NLongi==0) rownames(x$SurvEff[[1]]) <- sapply(rownames(x$SurvEff[[1]]), function(x) gsub("_S1", "", x))
        cat("\nSurvival outcome\n")
      }else{
        cat(paste0("\nSurvival outcome (S", i, ")\n"))
      }
      print(round(x$SurvEff[[i]], 4))
      if(!is.null(x$ReffListS[[i]])){
        if(x$NSurv==1){
          if(x$NLongi==0) rownames(x$ReffListS[[1]]) <- sapply(rownames(x$ReffListS[[1]]), function(x) gsub("_S1", "", x))
          if(!x$sdcor) cat(paste0("\nFrailty term variance\n")) else cat(paste0("\nFrailty term standard deviation\n"))
        }else{
          if(!x$sdcor) cat(paste0("\nFrailty term variance (S", i, ")\n")) else cat(paste0("\nFrailty term standard deviation (S", i, ")\n"))
        }
        print(round(x$ReffListS[[i]], 4))
      }
    }
  }
  if(!is.null(x$AssocLS)){
    if(dim(x$AssocLS)[1]>0){
      cat("\nAssociation longitudinal - survival\n")
      print(round(x$AssocLS, 4))
    }
  }
  if(!is.null(x$AssocSS)){
    if(dim(x$AssocSS)[1]>0){
      cat("\nAssociation survival - survival\n")
      print(round(x$AssocSS, 4))
    }
  }
  cat("\n")
  print(x$mlik[,1])
  if(!is.null(x$dic)) cat(paste0("\nDeviance Information Criterion: "), x$dic)
  if(!is.null(x$waic)) cat(paste0("\nWidely applicable Bayesian information criterion: "), x$waic)
  cat(paste0("\nComputation time: ", round(x$cpu.used[4], 2), " seconds"))
  cat("\n")
  invisible(x)
}
