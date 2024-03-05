#' @export

print.INLAjoint <- function(x, ...){
  cat(paste0("The model includes ", length(x$longOutcome), " longitudinal outcome(s) and ",
             length(x$survOutcome), " survival outcome(s). \n"))
  Nlongi <- length(x$longOutcome)
  if(Nlongi>0){
    for(i in 1:Nlongi){
      cat(paste0("Longitudinal outcome ", i, " (L", i, "):\n ", x$longOutcome[i], " (", x$famLongi[i], ")\n"))
    }
  }
  Nsurv <- length(x$survOutcome)
  if(Nsurv>0){
    for(i in 1:Nsurv){
      cat(paste0("Survival outcome ", i, " (S", i, "):\n Event time = ", x$SurvInfo[[i]]$nameTimeSurv,
                 "\n Event indicator = ", x$SurvInfo[[i]]$survOutcome, "\n baseline risk = ", x$basRisk[[i]], "\n"))
    }
  }
  Nassoc <- length(x$assoc)
  if(Nassoc>0){
    cat(paste0("There are  ", Nassoc, " association terms: \n"))
    for(i in 1:Nassoc){
      if(length(grep("NL_CV_CS", x$assoc[i]))==1){
        cat(paste0(" Non-linear current value and current slope: ", x$assoc[i], "\n"))
      }else if(length(grep("NL_CV", x$assoc[i]))==1){
        cat(paste0(" Non-linear current value: ", x$assoc[i], "\n"))
      }else if(length(grep("NL_CS", x$assoc[i]))==1){
        cat(paste0(" Non-linear current slope: ", x$assoc[i], "\n"))
      }else if(length(grep("NL_SRE", x$assoc[i]))==1){
        cat(paste0(" Non-linear shared random effects: ", x$assoc[i], "\n"))
      }else if(length(grep("CV_CS", x$assoc[i]))==1){
        cat(paste0(" Current value and current slope: ", x$assoc[i], "\n"))
      }else if(length(grep("CV", x$assoc[i]))==1){
        cat(paste0(" Current value: ", x$assoc[i], "\n"))
      }else if(length(grep("CS", x$assoc[i]))==1){
        cat(paste0(" Current slope: ", x$assoc[i], "\n"))
      }else if(length(grep("SRE_ind", x$assoc[i]))==1){
        cat(paste0(" Shared random effects (time fixed, each scaled independently): ", x$assoc[i], "\n"))
      }else if(length(grep("SRE", x$assoc[i]))==1){
        cat(paste0(" Shared random effects (time-dependent individual deviation): ", x$assoc[i], "\n"))
      }
    }
  }
}
