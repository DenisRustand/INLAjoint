#' @export



summary.INLAjoint <- function(obj, sdcor=FALSE, ...){
  if (!class(obj)=="INLAjoint"){
    stop("Please provide an object of class 'INLAjoint' (obtained with joint() function).\n")
  }
  out <- NULL
  class(obj) <- "inla"

  m.lstat.1 <- function(m) { #SD
    m <- inla.smarginal(m)
    ab <- inla.qmarginal(c(0.001, 0.999), m)
    ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
    m$x <- m$x[ii]
    m$y <- m$y[ii]
    moments <- inla.emarginal(function(lx) c(exp(-lx/2), exp(-lx)), m)
    q = exp(-inla.qmarginal(c(0.025, 0.5, 0.975), m)/2)
    return(list(mean = moments[1], sd = sqrt(max(0, moments[2]-moments[1]^2)), "0.025quant"=q[3], "0.5quant"=q[2], "0.975quant"=q[1]))
  }
  m.lstat.2 <- function(m) { #Variance
    m <- inla.smarginal(m)
    ab <- inla.qmarginal(c(0.001, 0.999), m)
    ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
    m$x <- m$x[ii]
    m$y <- m$y[ii]
    moments <- inla.emarginal(function(lx) c(exp(-lx), exp(-2*lx)), m)
    q = exp(-inla.qmarginal(c(0.025, 0.5, 0.975), m))
    return(list(mean = moments[1], sd = sqrt(max(0, moments[2]-moments[1]^2)), "0.025quant"=q[3], "0.5quant"=q[2], "0.975quant"=q[1]))
  }

  CompoFixed <- substring(obj$names.fixed, nchar(obj$names.fixed)-1, nchar(obj$names.fixed))
  Ncompo <- length(unique(CompoFixed))
  Mark <- unique(CompoFixed)
  Lmark <- Mark[which(substring(Mark, nchar(Mark)-1, nchar(Mark)-1)=="L")] # Longitudinal marker(s)
  Smark <- Mark[which(substring(Mark, nchar(Mark)-1, nchar(Mark)-1)=="S")] # Survival outcome(s)
  NLongi <- length(unique(substring(Lmark, nchar(Lmark), nchar(Lmark))))
  NSurv <- length(unique(substring(Smark, nchar(Smark), nchar(Smark))))
  BH_temp <- obj$summary.hyperpar[which(substring(rownames(obj$summary.hyperpar), nchar(rownames(obj$summary.hyperpar))-5, nchar(rownames(obj$summary.hyperpar)))=="hazard"), -which(colnames(obj$summary.hyperpar)=="mode")]
  if(!is.null(dim(BH_temp)[1])) BH_temp2 <- vector("list", dim(BH_temp)[1])
  BH <- NULL
  if(dim(BH_temp)[1]==1){
    if(!sdcor){
      BH_temp2[[1]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
      }, error = function(error_message) {
        message("Warning: there is an issue with baseline risk variance, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
        message(error_message)
        BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
      })
      BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
      rownames(BH) <- "Baseline (var)"
    }else{
      BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
      }, error = function(error_message) {
        message("Warning: there is an issue with baseline risk standard deviation, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
        message(error_message)
        BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
      })
      BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
      rownames(BH) <- "Baseline (sd)"
    }
  }else if(dim(BH_temp)[1]>1){
    for(i in 1:dim(BH_temp)[1]){
      if(!sdcor){
        BH_temp2[[1]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp)[i], 2, nchar(rownames(BH_temp)[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with baseline risk variance, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
          message(error_message)
          BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
        BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
        rownames(BH)[i] <- paste0("Baseline_S", i, " (var)")
      }else{
        BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp)[i], 2, nchar(rownames(BH_temp)[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with baseline risk standard deviation, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
          message(error_message)
          BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
        BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
        rownames(BH)[i] <- paste0("Baseline_S", i, " (sd)")
      }
    }
  }
  if(!is.null(BH)) colnames(BH) <- c("mean", "sd", "0.025quant", "0.5quant","0.975quant")
  ResErrNames <- rownames(obj$summary.hyperpar)[which(substring(rownames(obj$summary.hyperpar), 1, 26)%in%c("Precision for the Gaussian", "Precision for the lognorma"))]
  VarErr <- vector("list", length(ResErrNames))
  if(length(ResErrNames)>0){
    for(i in 1:length(ResErrNames)){
      if(!sdcor){
        VarErr[[i]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(ResErrNames[i], 2, nchar(ResErrNames[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with the variance of the residual error, you can try rerunning or scale the marker to fix it. It has been set to zero for now.\n")
          message(error_message)
          VarErr[[i]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
      }else{
        VarErr[[i]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(ResErrNames[i], 2, nchar(ResErrNames[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with the standard deviation of the residual error, you can try rerunning or scale the marker to fix it. It has been set to zero for now.\n")
          message(error_message)
          VarErr[[i]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
      }
    }
  }
  RandEff <- obj$summary.hyperpar[which(substring(rownames(obj$summary.hyperpar), 1, 5)=="Theta" |
                                          substring(rownames(obj$summary.hyperpar), 1, 16)=="Precision for ID"),]
  NRand <- length(unique(substring(rownames(RandEff), nchar(rownames(RandEff))-1, nchar(rownames(RandEff)))))
  AssocLS <- obj$summary.hyperpar[which(substring(rownames(obj$summary.hyperpar), 1, 4)=="Beta"), -which(colnames(obj$summary.hyperpar)=="mode")]
  if(dim(AssocLS)[1]>0) rownames(AssocLS) <- sapply(strsplit(rownames(AssocLS), "Beta for "), function(x) x[2])
  out$AssocLS <- AssocLS

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  if(NRand>0){
    NREcur <- 1
    ReffList <- vector("list", NRand)
    for(i in 1:NRand){
      RandEffi <- RandEff[which(substring(rownames(RandEff), nchar(rownames(RandEff)), nchar(rownames(RandEff)))==i),]
      NRandEffi <- dim(RandEffi)[1]
      NameRandEffi <- strsplit(rownames(RandEffi)[1], "for ")[[1]][2]
      if(NRandEffi==1){
        if(!sdcor){
          if(TRUE %in% c(c("Inf", "NaN") %in% RandEffi)){ # in case of infinite or not a number in the random effect hyperparameter
            Varmar <- RandEffi[1,]
          }else{
            Varmar <- m.lstat.2(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log precision for ", NameRandEffi, "`"))))
          }
        }else{
          if(TRUE %in% c(c("Inf", "NaN") %in% RandEffi)){ # in case of infinite or not a number in the random effect hyperparameter
            Varmar <- RandEffi[1,]
          }else{
            Varmar <- m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log precision for ", NameRandEffi, "`"))))
          }
        }
        ReffList[[i]] <- cbind("mean" = Varmar$mean,
                               "sd" = Varmar$sd,
                               "0.025quant" = Varmar$`0.025quant`,
                               "0.5quant" = Varmar$`0.5quant`,
                               "0.975quant" = Varmar$`0.975quant`)
        rownames(ReffList[[i]]) <- obj$REstruc[NREcur]
        NREcur <- NREcur + 1
      }else{
        NRE = (-1+sqrt(8*NRandEffi+1))/2 # get number of rancom effects from length of Cholesky terms
        if(!sdcor){
          MC_samples <- inla.iidkd.sample(10^4, obj, NameRandEffi, return.cov=TRUE)
        }else{
          MC_samples <- inla.iidkd.sample(10^4, obj, NameRandEffi, return.cov=FALSE)
        }
        VarCov <- matrix(unlist(MC_samples), nrow = NRE^2)
        VarCovMeans <- matrix(rowMeans(VarCov),NRE,NRE)
        VarCovSD <- matrix(apply(VarCov, 1, sd),NRE,NRE)
        VarCov2.5 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)),NRE,NRE)# 2.5% lower CI
        VarCov50 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)),NRE,NRE)# 50% quantile
        VarCov97.5 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)),NRE,NRE)# 97.5% upper CI
        #VarCovMode <- matrix(apply(VarCov, 1, Mode),NRE,NRE)

        ReffList[[i]] <- cbind("mean" = c(diag(VarCovMeans), VarCovMeans[lower.tri(VarCovMeans)]),
                               "sd" = c(diag(VarCovSD), VarCovSD[lower.tri(VarCovSD)]),
                               "0.025quant" = c(diag(VarCov2.5), VarCov2.5[lower.tri(VarCov2.5)]),
                               "0.5quant" = c(diag(VarCov50), VarCov50[lower.tri(VarCov50)]),
                               "0.975quant" = c(diag(VarCov97.5), VarCov97.5[lower.tri(VarCov97.5)]))#,
        #"mode" = round(c(diag(VarCovMode), VarCovMode[lower.tri(VarCovMode)]),9))
        NamesRE <- obj$REstruc[NREcur:(NREcur+(NRE-1))]
        triRE <- NULL # set names for covariance parameters
        for(j in 1:(length(NamesRE)-1)){
          for(k in (j+1):(length(NamesRE))){
            triRE <- c(triRE, paste0(NamesRE[j], ":", NamesRE[k]))
          }
        }
        NREcur <- NREcur + NRE
        rownames(ReffList[[i]]) <- c(NamesRE, triRE)
      }
    }
    out$ReffList <- ReffList
  }
  if(!is.null(obj$famLongi)){
    Nerr <- 1 #  identify error terms
    out$famLongi <- obj$famLongi
    FixedEff <- vector("list", NLongi)
    for(i in 1:NLongi){
      FixedEffi <- obj$summary.fixed[which(substring(rownames(obj$summary.fixed), nchar(rownames(obj$summary.fixed))-1, nchar(rownames(obj$summary.fixed)))==paste0("L", i)), -which(colnames(obj$summary.fixed)%in%c("mode","kld"))]
      rownames(FixedEffi) <- gsub("\\.X\\.", ":", rownames(FixedEffi))
      if(obj$famLongi[i] %in% c("gaussian", "lognormal")){
        if(!sdcor){
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (var)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }else{
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (sd)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }
        Nerr <- Nerr+1
      }else{
        FixedEff[[i]] <- FixedEffi
      }
    }
    out$FixedEff <- FixedEff
  }

  if(NSurv>0){
    BH <- as.data.frame(BH)
    SurvEff <- vector("list", NSurv)
    for(i in 1:NSurv){
      SurvEffi <- rbind(BH[i,], obj$summary.fixed[which(substring(rownames(obj$summary.fixed), nchar(rownames(obj$summary.fixed))-1, nchar(rownames(obj$summary.fixed)))==paste0("S", i)), -which(colnames(obj$summary.fixed)%in%c("mode","kld"))])
      rownames(SurvEffi) <- gsub("\\.X\\.", ":", rownames(SurvEffi))
      SurvEff[[i]] <- SurvEffi
    }
    out$SurvEff <- SurvEff
  }
  out$sdcor <- sdcor
  out$dic <- obj$dic$dic
  out$waic <- obj$waic$waic
  out$cpo <- obj$cpo$cpo
  out$gcpo <- obj$gcpo$gcpo
  out$pit <- obj$cpo$pit
  out$NLongi <- NLongi
  out$NSurv <- NSurv
  out$NRand <- NRand
  out$mlik <- obj$mlik
  out$cpu.used <- obj$cpu.used
  class(out) <- "summary.INLAjoint"
  out
}
