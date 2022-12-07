#' @export



summary.INLAjoint <- function(obj, sdcor=FALSE, hazr=FALSE, ...){
  if (!"INLAjoint" %in% class(obj)){
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
  m.lstat.3 <- function(m) { #frailty SD
    m <- inla.smarginal(m)
    ab <- inla.qmarginal(c(0.001, 0.999), m)
    ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
    m$x <- m$x[ii]
    m$y <- m$y[ii]
    moments <- inla.emarginal(function(lx) c(exp(-lx/2), exp(-lx)), m)
    q = exp(-inla.qmarginal(c(0.025, 0.5, 0.975), m))
    return(list(mean = moments[1], sd = sqrt(max(0, moments[2]-moments[1]^2)), "0.025quant"=q[3], "0.5quant"=q[2], "0.975quant"=q[1]))
  }
  m.lstat.4 <- function(m) { #frailty variance
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
  BH <- NULL
  if(!is.null(dim(BH_temp)[1])){
    BH_temp2 <- vector("list", dim(BH_temp)[1])
    if(dim(BH_temp)[1]==1){
      if(!sdcor){
        BH_temp2[[1]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with baseline risk variance, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
          message(error_message)
          BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
        BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
        # rownames(BH) <- "Baseline risk (variance)_S1"
      }else{
        BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with baseline risk standard deviation, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
          message(error_message)
          BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
        BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
        # rownames(BH) <- "Baseline risk (sd)_S1"
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
          # rownames(BH)[i] <- paste0("Baseline risk (variance)_S", i)
        }else{
          BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp)[i], 2, nchar(rownames(BH_temp)[i])), "`"))))
          }, error = function(error_message) {
            message("Warning: there is an issue with baseline risk standard deviation, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
            message(error_message)
            BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
          })
          BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
          # rownames(BH)[i] <- paste0("Baseline risk (sd)_S", i)
        }
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
  REidentify <- which(substring(rownames(obj$summary.hyperpar), 1, 5)=="Theta" | # extract random effects from hyperparameters
                        (substring(rownames(obj$summary.hyperpar), 1, 16)=="Precision for ID"))
  REidentifyL <- REidentify[!REidentify %in% grep("_S",rownames(obj$summary.hyperpar))] # long
  REidentifyS <- REidentify[REidentify %in% grep("_S",rownames(obj$summary.hyperpar))] # surv
  RandEff <- obj$summary.hyperpar[REidentifyL,]
  RandEffS <- obj$summary.hyperpar[REidentifyS,]
  NRand <- length(unique(substring(rownames(RandEff), nchar(rownames(RandEff))-1, nchar(rownames(RandEff)))))
  NRandS <- length(unique(substring(rownames(RandEffS), nchar(rownames(RandEffS))-1, nchar(rownames(RandEffS)))))
  AssocALL <- obj$summary.hyperpar[which(substring(rownames(obj$summary.hyperpar), 1, 4)=="Beta"), -which(colnames(obj$summary.hyperpar)=="mode")]
  if(!is.null(AssocALL)){
    if(dim(AssocALL)[1]>0){
      if(hazr){
        for(j in 1:dim(AssocALL)[1]){ # hazards ratios
          m <- inla.smarginal(obj$marginals.hyperpar[[rownames(AssocALL)[j]]])
          ab <- inla.qmarginal(c(0.001, 0.999), m)
          ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
          m$x <- m$x[ii]
          m$y <- m$y[ii]
          trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)
          AssocALL[j, "mean"] <- trsf$mean
          AssocALL[j, "sd"] <- trsf$sd
          AssocALL[j, "0.025quant"] <- trsf$quant0.025
          AssocALL[j, "0.5quant"] <- trsf$quant0.5
          AssocALL[j, "0.975quant"] <- trsf$quant0.975
        }
        colnames(AssocALL)[1] <- "exp(mean)"
      }
      if(dim(AssocALL)[1]>0) rownames(AssocALL) <- sapply(strsplit(rownames(AssocALL), "Beta for "), function(x) x[2])
    }
    AssocLS <- AssocALL[which(substring(rownames(AssocALL), nchar(rownames(AssocALL))-5, nchar(rownames(AssocALL))-4)=="_L"),]
    #rownames(AssocLS) <- sapply(strsplit(rownames(AssocLS), "ID"), function(x) x[2])
    AssocSS <-AssocALL[which(substring(rownames(AssocALL), nchar(rownames(AssocALL))-5, nchar(rownames(AssocALL))-4)=="_S"),]
    #rownames(AssocSS) <- sapply(strsplit(rownames(AssocSS), "ID"), function(x) x[2])
  }
  out$AssocLS <- AssocLS
  out$AssocSS <- AssocSS
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
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (variance)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }else{
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (sd)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }
        Nerr <- Nerr+1
      }else if (obj$famLongi[i]=="pom"){
        FixedEff[[i]] <- rbind(obj$summary.hyperpar[grep("POM", rownames(obj$summary.hyperpar)),-6], FixedEffi)
      }else if(obj$famLongi[i]=="nbinomial"){
        FixedEff[[i]] <- rbind(obj$summary.hyperpar[grep("nbinomial", rownames(obj$summary.hyperpar)),-6], FixedEffi)
      }else if(obj$famLongi[i]=="Betabinomial"){
        FixedEff[[i]] <- rbind(obj$summary.hyperpar[grep("betabinomial", rownames(obj$summary.hyperpar)),-6], FixedEffi)
      }else if(obj$famLongi[i]=="gpoisson"){
        FixedEff[[i]] <- rbind(obj$summary.hyperpar[grep("gpoisson", rownames(obj$summary.hyperpar)),-6], FixedEffi)
      }else if(length(grep("zeroinflated", obj$famLongi[i]))>0){
        FixedEff[[i]] <- rbind(obj$summary.hyperpar[grep("zero-inflated", rownames(obj$summary.hyperpar)),-6], FixedEffi)
      }else{
        FixedEff[[i]] <- FixedEffi
      }
    }
    out$FixedEff <- FixedEff
  }
  if(NSurv>0){
    BH <- as.data.frame(BH)
    BHW <- vector("list", NSurv) # baseline risk
    SurvEff <- vector("list", NSurv)
    BaselineValues <- vector("list", NSurv) # values of piecewise constant baseline
    nbl <- 1 # to keep track of baseline risk in case of multiple survival outcomes
    nbl2 <- 1 # to keep track of baseline risk in case of multiple survival outcomes
    MCure <- NULL # Mixture cure part
    for(i in 1:NSurv){
      if(obj$basRisk[[i]]=="exponentialsurv"){
        nameRisk <- "Exponential (rate)"
      }else if(obj$basRisk[[i]]=="weibullsurv"){
        nameRisk <- "Weibull (scale)"
        BHW[[i]] <- obj$summary.hyperpar[grep("weibull", rownames(obj$summary.hyperpar)),][nbl2,-6]
        rownames(BHW[[i]]) <- paste0("Weibull (shape)_S", i)
        nbl2 <- nbl2+1
      }else{
        nameRisk <- "Baseline risk (mean)"
        BHW[[i]] <- BH[nbl,]
        rownames(BHW[[i]]) <- paste0("Baseline risk (variance)_S",i)
        nbl <- nbl+1
        BHmean <- NULL
        BSsd <- NULL
        BHlo <- NULL
        BHme <- NULL
        BHup <- NULL
        for(j in 1:length(obj$marginals.random[[i]])){
          # m <- inla.smarginal(obj$marginals.random[[i]][[j]])
          # ab <- inla.qmarginal(c(0.001, 0.999), m)
          # ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
          # m$x <- m$x[ii]
          # m$y <- m$y[ii]
          # m <- inla.smarginal(m)
          #trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)
          #BHmean <- c(BHmean, trsf$mean)
          #BSsd <- c(BSsd, trsf$sd)
          Mm <- inla.qmarginal(c(0.025, 0.5, 0.975), obj$marginals.random[[i]][[j]])
          BHlo <- c(BHlo, exp(Mm[1]))
          BHme <- c(BHme, exp(Mm[2]))
          BHup <- c(BHup, exp(Mm[3]))
        }
        BaselineValues[[i]] <- cbind(time=obj$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                                     # mean=BHmean,
                                     # sd=BSsd,
                                     lower=BHlo,
                                     median=BHme,
                                     upper=BHup)
      }
      if(obj$basRisk[[i]] %in% c("exponentialsurv", "weibullsurv") & !is.null(obj$cureVar[[i]])){
        MCure <- obj$summary.hyperpar[grep("Weibull-Cure", rownames(obj$summary.hyperpar)),-6]
        rownames(MCure) <- obj$cureVar[[i]]
      }
      SurvEffi <- rbind(BHW[[i]], obj$summary.fixed[which(substring(rownames(obj$summary.fixed), nchar(rownames(obj$summary.fixed))-1, nchar(rownames(obj$summary.fixed)))==paste0("S", i)), -which(colnames(obj$summary.fixed)%in%c("mode","kld"))])
      rownames(SurvEffi) <- gsub("\\.X\\.", ":", rownames(SurvEffi))
      rownames(SurvEffi)[grep("Intercept", rownames(SurvEffi))] <- paste0(nameRisk, "_S", i)
      if(!is.null(obj$marginals.fixed[[paste0("Intercept_S",i)]])){
        m <- inla.smarginal(obj$marginals.fixed[[paste0("Intercept_S",i)]]) # baseline risk mean
        ab <- inla.qmarginal(c(0.001, 0.999), m)
        ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
        m$x <- m$x[ii]
        m$y <- m$y[ii]
        trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)

        SurvEffi[paste0(nameRisk, "_S", i), "mean"] <- trsf$mean
        SurvEffi[paste0(nameRisk, "_S", i), "sd"] <- trsf$sd
        SurvEffi[paste0(nameRisk, "_S", i), "0.025quant"] <- trsf$quant0.025
        SurvEffi[paste0(nameRisk, "_S", i), "0.5quant"] <- trsf$quant0.5
        SurvEffi[paste0(nameRisk, "_S", i), "0.975quant"] <- trsf$quant0.975
      }
      if(hazr){
        for(j in 1:dim(SurvEffi)[1]){ # hazards ratios
          if(!j%in%c(grep("aseline", rownames(SurvEffi)),
                     grep("(rate)", rownames(SurvEffi)),
                     grep("(shape)", rownames(SurvEffi)),
                     grep("(scale)", rownames(SurvEffi)))){
            RNM <- gsub(":", "\\.X\\.", rownames(SurvEffi)[j])
            m <- inla.smarginal(obj$marginals.fixed[[RNM]])
            ab <- inla.qmarginal(c(0.001, 0.999), m)
            ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
            m$x <- m$x[ii]
            m$y <- m$y[ii]
            trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)
            SurvEffi[j, "mean"] <- trsf$mean
            SurvEffi[j, "sd"] <- trsf$sd
            SurvEffi[j, "0.025quant"] <- trsf$quant0.025
            SurvEffi[j, "0.5quant"] <- trsf$quant0.5
            SurvEffi[j, "0.975quant"] <- trsf$quant0.975
          }
        }
      }
      if(!is.null(MCure)){
        if(hazr) rownames(MCure) <- paste0(rownames(MCure), " (not exp!)")
        SurvEffi <- rbind(MCure, SurvEffi)
        if(hazr) colnames(SurvEffi)[1] <- "exp(mean)"
      }else if(hazr){
        colnames(SurvEffi)[1] <- "exp(mean)"
      }
      SurvEff[[i]] <- SurvEffi
    }
    out$SurvEff <- SurvEff

    if(NRandS>0){
      NREcurS <- 1
      ReffListS <- vector("list", NSurv)
      for(i in 1:NRandS){
        RandEffiS <- RandEffS[which(substring(rownames(RandEffS), nchar(rownames(RandEffS)), nchar(rownames(RandEffS)))==i),]
        NRandEffiS <- dim(RandEffiS)[1]
        NameRandEffiS <- strsplit(rownames(RandEffiS)[1], "for ")[[1]][2]
        if(NRandEffiS==1){
          if(!sdcor){
            if(TRUE %in% c(c("Inf", "NaN") %in% RandEffiS)){ # in case of infinite or not a number in the random effect hyperparameter
              VarmarS <- RandEffiS[1,]
            }else{
              VarmarS <- m.lstat.4(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log precision for ", NameRandEffiS, "`"))))
            }
          }else{
            if(TRUE %in% c(c("Inf", "NaN") %in% RandEffiS)){ # in case of infinite or not a number in the random effect hyperparameter
              VarmarS <- RandEffiS[1,]
            }else{
              VarmarS <- m.lstat.3(eval(parse(text=paste0("obj$internal.marginals.hyperpar$`Log precision for ", NameRandEffiS, "`"))))
            }
          }
          ReffListS[[i]] <- cbind("mean" = VarmarS$mean,
                                 "sd" = VarmarS$sd,
                                 "0.025quant" = VarmarS$`0.025quant`,
                                 "0.5quant" = VarmarS$`0.5quant`,
                                 "0.975quant" = VarmarS$`0.975quant`)
          rownames(ReffListS[[i]]) <- obj$REstrucS[NREcurS]
          NREcurS <- NREcurS + 1
        }
      }
      out$ReffListS <- ReffListS
    }
  }
  if(NSurv>0) out$BaselineValues <- BaselineValues
  out$sdcor <- sdcor
  out$dic <- obj$dic$dic
  out$waic <- obj$waic$waic
  out$cpo <- obj$cpo$cpo
  out$gcpo <- obj$gcpo$gcpo
  out$pit <- obj$cpo$pit
  out$NLongi <- NLongi
  out$NSurv <- NSurv
  out$NRand <- NRand
  out$NRandS <- NRandS
  out$mlik <- obj$mlik
  out$cpu.used <- obj$cpu.used
  class(out) <- "summary.INLAjoint"
  out
}
