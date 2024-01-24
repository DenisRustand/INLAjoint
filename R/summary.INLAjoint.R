#' @export



summary.INLAjoint <- function(object, ...){
  arguments <- list(...)
  if(!is.null(arguments$hazr)) hr=arguments$hazr
  if(is.null(arguments$hr)) hr=F else hr=arguments$hr
  if(is.null(arguments$sdcor)) sdcor=F else sdcor=arguments$sdcor
  if(is.null(arguments$NsampRE)) NsampRE=2000 else NsampRE=arguments$NsampRE
  if (!"INLAjoint" %in% class(object)){
    stop("Please provide an object of class 'INLAjoint' (obtained with joint() function).\n")
  }
  out <- NULL
  class(object) <- "inla"
  m.lstat.1 <- function(m) { #SD
    m <- INLA::inla.smarginal(m)
    ab <- INLA::inla.qmarginal(c(0.001, 0.999), m)
    ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
    m$x <- m$x[ii]
    m$y <- m$y[ii]
    moments <- INLA::inla.emarginal(function(lx) c(exp(-lx/2), exp(-lx)), m)
    q = exp(-INLA::inla.qmarginal(c(0.025, 0.5, 0.975), m)/2)
    return(list(mean = moments[1], sd = sqrt(max(0, moments[2]-moments[1]^2)), "0.025quant"=q[3], "0.5quant"=q[2], "0.975quant"=q[1]))
  }
  m.lstat.2 <- function(m) { #Variance
    m <- INLA::inla.smarginal(m)
    ab <- INLA::inla.qmarginal(c(0.001, 0.999), m)
    ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
    m$x <- m$x[ii]
    m$y <- m$y[ii]
    moments <- INLA::inla.emarginal(function(lx) c(exp(-lx), exp(-2*lx)), m)
    q = exp(-INLA::inla.qmarginal(c(0.025, 0.5, 0.975), m))
    return(list(mean = moments[1], sd = sqrt(max(0, moments[2]-moments[1]^2)), "0.025quant"=q[3], "0.5quant"=q[2], "0.975quant"=q[1]))
  }
  CompoFixed <- gsub("_","", substring(object$names.fixed, nchar(object$names.fixed)-2, nchar(object$names.fixed)))
  Ncompo <- length(unique(CompoFixed))
  Mark <- unique(CompoFixed)
  Lmark <- Mark[grep("L", Mark)] # Longitudinal marker(s)
  Smark <- Mark[grep("S", Mark)] # Survival outcome(s)
  NLongi <- length(unique(Lmark))
  NSurv <- length(unique(Smark))
  Hnames <- rownames(object$summary.hyperpar)
  BH_temp <- object$summary.hyperpar[which(substring(Hnames, nchar(Hnames)-5, nchar(Hnames))=="hazard"), -which(colnames(object$summary.hyperpar)=="mode")]
  BH <- NULL
  if(!is.null(dim(BH_temp)[1])){
    BH_temp2 <- vector("list", dim(BH_temp)[1])
    if(dim(BH_temp)[1]==1){
      if(!sdcor){
        BH_temp2[[1]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with baseline risk variance, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
          message(error_message)
          BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
        BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
        # rownames(BH) <- "Baseline risk (variance)_S1"
      }else{
        BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp), 2, nchar(rownames(BH_temp))), "`"))))
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
          BH_temp2[[1]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp)[i], 2, nchar(rownames(BH_temp)[i])), "`"))))
          }, error = function(error_message) {
            message("Warning: there is an issue with baseline risk variance, you can try rerunning, scale the event times or change the number of intervals in the baseline risk to fix it. It has been set to zero for now.\n")
            message(error_message)
            BH_temp2[[1]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
          })
          BH <- rbind(BH, unlist(c(BH_temp2[[1]]["mean"], BH_temp2[[1]]["sd"], BH_temp2[[1]]["0.025quant"], BH_temp2[[1]]["0.5quant"], BH_temp2[[1]]["0.975quant"])))
          # rownames(BH)[i] <- paste0("Baseline risk (variance)_S", i)
        }else{
          BH_temp2[[1]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(rownames(BH_temp)[i], 2, nchar(rownames(BH_temp)[i])), "`"))))
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
  ResErrNames <- Hnames[which(substring(Hnames, 1, 26)%in%c("Precision for the Gaussian", "Precision for the lognorma"))]
  VarErr <- vector("list", length(ResErrNames))
  if(length(ResErrNames)>0){
    for(i in 1:length(ResErrNames)){
      if(!sdcor){
        VarErr[[i]] <- tryCatch({m.lstat.2(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(ResErrNames[i], 2, nchar(ResErrNames[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with the variance of the residual error, you can try rerunning or scale the marker to fix it. It has been set to zero for now.\n")
          message(error_message)
          VarErr[[i]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
      }else{
        VarErr[[i]] <- tryCatch({m.lstat.1(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log p", substring(ResErrNames[i], 2, nchar(ResErrNames[i])), "`"))))
        }, error = function(error_message) {
          message("Warning: there is an issue with the standard deviation of the residual error, you can try rerunning or scale the marker to fix it. It has been set to zero for now.\n")
          message(error_message)
          VarErr[[i]] <- c("mean"=0, "sd"=0, "0.025quant"=0, "0.5quant"=0, "0.975quant"=0)
        })
      }
    }
  }
  REidentify <- which(substring(Hnames, 1, 5)=="Theta" | # extract random effects from hyperparameters
                        (substring(Hnames, 1, 16)=="Precision for ID"))
  REidentifyL <- REidentify[!REidentify %in% grep("_S",Hnames)] # long
  REidentifyS <- REidentify[REidentify %in% grep("_S",Hnames)] # surv
  RandEff <- object$summary.hyperpar[REidentifyL,]
  RandEffS <- object$summary.hyperpar[REidentifyS,]
  if(!is.null(object$corLong)){
    NRand <- ifelse(object$corLong, 1, length(object$longOutcome))#length(unique(substring(rownames(RandEff), nchar(rownames(RandEff))-1, nchar(rownames(RandEff)))))
  }else{
    NRand <- 0
  }
  NRandS <- length(unique(substring(rownames(RandEffS), nchar(rownames(RandEffS))-1, nchar(rownames(RandEffS)))))
  AssocALL <- object$summary.hyperpar[which(substring(Hnames, 1, 5)=="Beta "), -which(colnames(object$summary.hyperpar)=="mode")]
  AssocNL <- object$summary.hyperpar[sort(c(grep("(scopy theta)", Hnames), grep("(scopy slope)", Hnames))), -which(colnames(object$summary.hyperpar)=="mode")]
  AssocLS <- NULL
  AssocSS <- NULL
  if(!is.null(AssocALL)){
    if(dim(AssocALL)[1]>0){
      if(hr){
        for(j in 1:dim(AssocALL)[1]){ # hazards ratios
          m <- INLA::inla.smarginal(object$marginals.hyperpar[[rownames(AssocALL)[j]]])
          ab <- INLA::inla.qmarginal(c(0.001, 0.999), m)
          ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
          m$x <- m$x[ii]
          m$y <- m$y[ii]
          trsf <- INLA::inla.zmarginal(INLA::inla.tmarginal(function(x) exp(x), m), silent=T)
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
    AssocLS <- AssocALL[intersect(grep("_L", rownames(AssocALL)), grep("_S", rownames(AssocALL))),]
    # AssocLS <- AssocALL[which(substring(rownames(AssocALL), nchar(rownames(AssocALL))-5, nchar(rownames(AssocALL))-4)=="_L"),]
    #rownames(AssocLS) <- sapply(strsplit(rownames(AssocLS), "ID"), function(x) x[2])
    AssocSS <-AssocALL[setdiff(grep("_S", rownames(AssocALL)), grep("_L", rownames(AssocALL))),]
    #rownames(AssocSS) <- sapply(strsplit(rownames(AssocSS), "ID"), function(x) x[2])
  }
  if(!is.null(AssocNL)){
    if(dim(AssocNL)[1]>0){
      # if(hr){
      #   for(j in 1:dim(AssocNL)[1]){ # hazards ratios
      #     m <- inla.smarginal(object$marginals.hyperpar[[rownames(AssocNL)[j]]])
      #     ab <- inla.qmarginal(c(0.001, 0.999), m)
      #     ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
      #     m$x <- m$x[ii]
      #     m$y <- m$y[ii]
      #     trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)
      #     AssocNL[j, "mean"] <- trsf$mean
      #     AssocNL[j, "sd"] <- trsf$sd
      #     AssocNL[j, "0.025quant"] <- trsf$quant0.025
      #     AssocNL[j, "0.5quant"] <- trsf$quant0.5
      #     AssocNL[j, "0.975quant"] <- trsf$quant0.975
      #   }
      #   colnames(AssocNL)[1] <- "exp(mean)"
      # }
      if(dim(AssocNL)[1]>0) rownames(AssocNL) <- sapply(rownames(AssocNL), function(x) strsplit(x, " \\(scopy")[[1]][1])
    }
    AssocLS <- rbind(AssocLS, AssocNL)
    #rownames(AssocLS) <- sapply(strsplit(rownames(AssocLS), "ID"), function(x) x[2])
    # AssocSS <-AssocNL[which(substring(rownames(AssocNL), nchar(rownames(AssocNL))-5, nchar(rownames(AssocNL))-4)=="_S"),]
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
      if(!is.null(object$fixRE[[i]])){
        if(object$fixRE[[i]]) fix_i <- TRUE else fix_i <- FALSE
      }else{
        fix_i <- FALSE
      }
      if(!fix_i){ #is.null(object$mat_k[[i]])
        if(i>9) shiftRE <- 3 else shiftRE <- 2
        RandEffi <- RandEff[which(substring(rownames(RandEff), nchar(rownames(RandEff))-shiftRE, nchar(rownames(RandEff)))==paste0("_L", i)),]
        NRandEffi <- dim(RandEffi)[1]
        NameRandEffi <- strsplit(rownames(RandEffi)[1], "for ")[[1]][2]
        if(NRandEffi==1 | !object$corRE[[i]]){
          # browser()
          for(j in 1:nrow(RandEffi)){
            NameRandEffi <- strsplit(rownames(RandEffi)[j], "for ")[[1]][2]
            if(!sdcor){
              if(TRUE %in% c(c("Inf", "NaN") %in% RandEffi)){ # in case of infinite or not a number in the random effect hyperparameter
                Varmar <- RandEffi[j,]
              }else{
                Varmar <- m.lstat.2(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log precision for ", NameRandEffi, "`"))))
              }
            }else{
              if(TRUE %in% c(c("Inf", "NaN") %in% RandEffi)){ # in case of infinite or not a number in the random effect hyperparameter
                Varmar <- RandEffi[j,]
              }else{
                Varmar <- m.lstat.1(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log precision for ", NameRandEffi, "`"))))
              }
            }
            Refi <- cbind("mean" = Varmar$mean,
                                   "sd" = Varmar$sd,
                                   "0.025quant" = Varmar$`0.025quant`,
                                   "0.5quant" = Varmar$`0.5quant`,
                                   "0.975quant" = Varmar$`0.975quant`)
            rownames(Refi) <- object$REstruc[NREcur]
            ReffList[[i]] <- rbind(ReffList[[i]], Refi)
            NREcur <- NREcur + 1
          }
        }else{
          NRE = (-1+sqrt(8*NRandEffi+1))/2 # get number of rancom effects from length of Cholesky terms
          if(!sdcor){
            MC_samples <- INLA::inla.iidkd.sample(NsampRE, object, NameRandEffi, return.cov=TRUE)
          }else{
            MC_samples <- INLA::inla.iidkd.sample(NsampRE, object, NameRandEffi, return.cov=FALSE)
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
          NamesRE <- object$REstruc[NREcur:(NREcur+(NRE-1))]
          triRE <- NULL # set names for covariance parameters
          for(j in 1:(length(NamesRE)-1)){
            for(k in (j+1):(length(NamesRE))){
              triRE <- c(triRE, paste0(NamesRE[j], ":", NamesRE[k]))
            }
          }
          NREcur <- NREcur + NRE
          rownames(ReffList[[i]]) <- c(NamesRE, triRE)
        }
      }else{
        NRE <- ifelse(is.null(dim(object$mat_k[[i]])[1]), 1, dim(object$mat_k[[i]])[1])
        if(sdcor){
          if(NRE>1){
            sdmat <- cov2cor(object$mat_k[[i]])
            diag(sdmat) <- sqrt(diag(object$mat_k[[i]]))
            ReffList[[i]] <- sdmat
          }else{
            ReffList[[i]] <- as.matrix(sqrt(object$mat_k[[i]]))
          }
        }else{
          ReffList[[i]] <- object$mat_k[[i]]
        }
        # if(NRE>1){
          rownames(ReffList[[i]]) <- colnames(ReffList[[i]]) <- object$REstruc[NREcur:(NREcur+(NRE-1))]
        # }else{
        #   names(ReffList[[i]]) <- object$REstruc[NREcur:(NREcur+(NRE-1))]
        # }
        NREcur <- NREcur + NRE
      }
    }
    out$ReffList <- ReffList
  }
  if(!is.null(object$famLongi) & NLongi>0){
    Nerr <- 1 #  identify error terms
    out$famLongi <- object$famLongi
    FixedEff <- vector("list", NLongi)
    for(i in 1:NLongi){
      if(i>9) shiftFE <- 2 else shiftFE <- 1
      FixedEffi <- object$summary.fixed[which(substring(rownames(object$summary.fixed), nchar(rownames(object$summary.fixed))-shiftFE, nchar(rownames(object$summary.fixed)))==paste0("L", i)), -which(colnames(object$summary.fixed)%in%c("mode","kld"))]
      rownames(FixedEffi) <- gsub("\\.X\\.", ":", rownames(FixedEffi))
      if(object$famLongi[i] %in% c("gaussian", "lognormal")){
        if(!sdcor){
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (variance)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }else{
          FixedEff[[i]] <- rbind(FixedEffi, "Res. err. (sd)" = c(VarErr[[Nerr]]["mean"], VarErr[[Nerr]]["sd"], VarErr[[Nerr]]["0.025quant"], VarErr[[Nerr]]["0.5quant"], VarErr[[Nerr]]["0.975quant"]))
        }
        Nerr <- Nerr+1
      }else if (object$famLongi[i]=="pom"){
        FixedEff[[i]] <- rbind(object$summary.hyperpar[grep("POM", Hnames),-6], FixedEffi)
      }else if(object$famLongi[i]=="nbinomial"){
        FixedEff[[i]] <- rbind(object$summary.hyperpar[grep("nbinomial", Hnames),-6], FixedEffi)
      }else if(object$famLongi[i]=="Betabinomial"){
        FixedEff[[i]] <- rbind(object$summary.hyperpar[grep("betabinomial", Hnames),-6], FixedEffi)
      }else if(object$famLongi[i]=="gpoisson"){
        FixedEff[[i]] <- rbind(object$summary.hyperpar[grep("gpoisson", Hnames),-6], FixedEffi)
      }else if(length(grep("zeroinflated", object$famLongi[i]))>0){
        FixedEff[[i]] <- rbind(object$summary.hyperpar[grep("zero-inflated", Hnames),-6], FixedEffi)
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
      if(object$basRisk[[i]]=="exponentialsurv"){
        nameRisk <- "Exponential (rate)"
      }else if(object$basRisk[[i]]=="weibullsurv"){
        nameRisk <- "Weibull (scale)"
        BHW[[i]] <- object$summary.hyperpar[grep("weibull", Hnames),][nbl2,-6]
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
        for(j in 1:length(object$marginals.random[[i]])){
          # m <- inla.smarginal(object$marginals.random[[i]][[j]])
          # ab <- inla.qmarginal(c(0.001, 0.999), m)
          # ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
          # m$x <- m$x[ii]
          # m$y <- m$y[ii]
          # m <- inla.smarginal(m)
          #trsf <- inla.zmarginal(inla.tmarginal(function(x) exp(x), m), silent=T)
          #BHmean <- c(BHmean, trsf$mean)
          #BSsd <- c(BSsd, trsf$sd)
          Mm <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), object$marginals.random[[i]][[j]])
          BHlo <- c(BHlo, exp(Mm[1]))
          BHme <- c(BHme, exp(Mm[2]))
          BHup <- c(BHup, exp(Mm[3]))
        }
        BaselineValues[[i]] <- cbind(time=object$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                                     # mean=BHmean,
                                     # sd=BSsd,
                                     lower=BHlo,
                                     median=BHme,
                                     upper=BHup)
      }
      if(object$basRisk[[i]] %in% c("exponentialsurv", "weibullsurv") & !is.null(object$cureVar[[i]])){
        MCure <- object$summary.hyperpar[grep("Weibull-Cure", Hnames),-6]
        rownames(MCure) <- object$cureVar[[i]]
      }
      SurvEffi <- rbind(BHW[[i]], object$summary.fixed[which(substring(rownames(object$summary.fixed), nchar(rownames(object$summary.fixed))-1, nchar(rownames(object$summary.fixed)))==paste0("S", i)), -which(colnames(object$summary.fixed)%in%c("mode","kld"))])
      rownames(SurvEffi) <- gsub("\\.X\\.", ":", rownames(SurvEffi))
      rownames(SurvEffi)[grep("Intercept", rownames(SurvEffi))] <- paste0(nameRisk, "_S", i)
      if(!is.null(object$marginals.fixed[[paste0("Intercept_S",i)]])){
        m <- INLA::inla.smarginal(object$marginals.fixed[[paste0("Intercept_S",i)]]) # baseline risk mean
        ab <- INLA::inla.qmarginal(c(0.001, 0.999), m)
        ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
        m$x <- m$x[ii]
        m$y <- m$y[ii]
        if(object$variant==0){
          trsf <- INLA::inla.zmarginal(INLA::inla.tmarginal(function(x) exp(x), m), silent=T)
        }else if(object$variant==1){
          #SMPLweib <- inla.posterior.sample(100, object)
          warning("See inla.doc('weibullsurv') for details about variant 1, the standard weibull baseline survival model uses variant=0")
          trsf <- INLA::inla.zmarginal(INLA::inla.tmarginal(function(x) exp(x), m), silent=T)

          #inla.posterior.sample.eval(function(...){Intercept_S1*alpha parameter for weibullsurv} , SMPLweib)
          # need to do the alpha*Beta0 for intercept ("scale" (should it be renamed?)) here and exp(alpha*betas) below for hazard ratios
          #}
        }

        SurvEffi[paste0(nameRisk, "_S", i), "mean"] <- trsf$mean
        SurvEffi[paste0(nameRisk, "_S", i), "sd"] <- trsf$sd
        SurvEffi[paste0(nameRisk, "_S", i), "0.025quant"] <- trsf$quant0.025
        SurvEffi[paste0(nameRisk, "_S", i), "0.5quant"] <- trsf$quant0.5
        SurvEffi[paste0(nameRisk, "_S", i), "0.975quant"] <- trsf$quant0.975
      }
      if(hr){
        for(j in 1:dim(SurvEffi)[1]){ # hazards ratios
          if(!j%in%c(grep("aseline", rownames(SurvEffi)),
                     grep("(rate)", rownames(SurvEffi)),
                     grep("(shape)", rownames(SurvEffi)),
                     grep("(scale)", rownames(SurvEffi)))){
            RNM <- gsub(":", "\\.X\\.", rownames(SurvEffi)[j])
            m <- INLA::inla.smarginal(object$marginals.fixed[[RNM]])
            ab <- INLA::inla.qmarginal(c(0.001, 0.999), m)
            ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
            m$x <- m$x[ii]
            m$y <- m$y[ii]
            trsf <- INLA::inla.zmarginal(INLA::inla.tmarginal(function(x) exp(x), m), silent=T)
            SurvEffi[j, "mean"] <- trsf$mean
            SurvEffi[j, "sd"] <- trsf$sd
            SurvEffi[j, "0.025quant"] <- trsf$quant0.025
            SurvEffi[j, "0.5quant"] <- trsf$quant0.5
            SurvEffi[j, "0.975quant"] <- trsf$quant0.975
            # ?inla.posterior.sample.eval
          }
        }
      }
      if(!is.null(MCure)){
        if(hr) rownames(MCure) <- paste0(rownames(MCure), " (not exp!)")
        SurvEffi <- rbind(MCure, SurvEffi)
        if(hr) colnames(SurvEffi)[1] <- "exp(mean)"
      }else if(hr){
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
              VarmarS <- m.lstat.2(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log precision for ", NameRandEffiS, "`"))))
            }
          }else{
            if(TRUE %in% c(c("Inf", "NaN") %in% RandEffiS)){ # in case of infinite or not a number in the random effect hyperparameter
              VarmarS <- RandEffiS[1,]
            }else{
              VarmarS <- m.lstat.1(eval(parse(text=paste0("object$internal.marginals.hyperpar$`Log precision for ", NameRandEffiS, "`"))))
            }
          }
          ReffListS[[i]] <- cbind("mean" = VarmarS$mean,
                                  "sd" = VarmarS$sd,
                                  "0.025quant" = VarmarS$`0.025quant`,
                                  "0.5quant" = VarmarS$`0.5quant`,
                                  "0.975quant" = VarmarS$`0.975quant`)
          rownames(ReffListS[[i]]) <- object$REstrucS[NREcurS]
          NREcurS <- NREcurS + 1
        }
      }
      out$ReffListS <- ReffListS
    }
  }
  if(NSurv>0) out$BaselineValues <- BaselineValues
  out$sdcor <- sdcor
  out$dic <- object$dic$dic
  out$waic <- object$waic$waic
  out$cpo <- object$cpo$cpo
  out$gcpo <- object$gcpo$gcpo
  out$pit <- object$cpo$pit
  out$NLongi <- NLongi
  out$NSurv <- NSurv
  out$NRand <- NRand
  out$NRandS <- NRandS
  out$mlik <- object$mlik
  out$cpu.used <- object$cpu.used
  class(out) <- "summary.INLAjoint"
  out
}
