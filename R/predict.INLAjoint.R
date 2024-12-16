#' Computes predictions for a given model fitted with INLAjoint
#'
#' @description This function allows to compute predictions for a given model fitted with INLAjoint,
#' the default behavior (without arguments) returns fitted values for each component of the model. It
#' is also possible to supply a dataset for which predictions are required, this dataset must have
#' the same structure as the dataset used for the model fitting (i.e., same columns). The default
#' returned predictions corresponds to the linear predictors for each outcomes.
#'
#' @param object an object that contains a model fitted with INLAjoint.
#' @param newData a dataset with the same columns as those used to fit the model. When using a longitudinal
#' marker to predict longitudinal and subsequent survival outcomes, only the longitudinal information (i.e.,
#' structure of the longitudinal data) is required. It is also possible to predict the average trajectories
#' conditional on covariates by setting the value of the longitudinal outcomes included in the model to NA.
#' @param newDataSurv a dataset for survival information (only useful when both longitudinal and survival
#' data are provided for the predictions, otherwise using the argument newData is working too).
#' @param timePoints a vector of the time points at which predictions are computed (for both longitudinal
#' and survival outcomes), this also control the precision of the integration for time-dependent shared
#' terms and the computation of cumulative risks (e.g., for survival or CIF curves), thus many time points
#' will increase the accuracy of predictions. Default is NULL as these time points are automatically computed
#' when not defined manually.
#' @param NtimePoints number of time points at which the predictions are computed (for both longitudinal
#' and survival outcomes), these time points are equidistant between time 0 and horizon time.
#' This also control the precision of the integration for time-dependent shared
#' terms and the computation of cumulative risks (e.g., for survival or CIF curves), thus many time points
#' will increase the accuracy of predictions.
#' @param NsampleHY number of samples for hyperparameters used to assess uncertainty
#' when computing predictions. Default is 20.
#' @param NsampleFE number of samples of fixed effects for each hyperparameters samples
#' used to assess uncertainty when computing predictions. Default is 30 (i.e., 30 x NsampleHY).
#' @param NsampleRE number of random effects realizations for each sample specified in 'NsampleHY' and
#' 'NsampleFE'. Default is 50 (i.e., 50 x NsampleFE x NsampleHY, resulting in 20000 random effects samples
#' per new individual with default values). These random effects realizations are conditional on
#' observed longitudinal outcomes values provided in 'newData' and survival time provided in 'newDataSurv'
#' when a survival model is included. If 'newDataSurv' is NULL, they are conditional on
#' survival up to latest longitudinal recorded measurement. When outcomes are
#' set to NA, the realizations are sampled from the marginal distribution of random effects.
#' @param id name of the individual id variable, default is NULL as it is automatically grabbed from the
#' fitted model but when fitting simple survival models, providing id when fitting the model is not
#' mandatory and thus this can be useful (an explicit message is printed in this specific case).
#' @param Csurv conditional survival, gives the starting value of the at-risk period (i.e., starting value
#' at which risk predictions for survival models are computed).
#' Default is the last longitudinal observation time provided in 'newData' but this is
#' replaced by the value of 'Csurv' when provided.
#' @param horizon horizon of the prediction.
#' @param baselineHaz method used to evaluate the baseline hazard value, default is 'interpolation'
#' which is currently recommended. Experimental alternatives are being developed, including 'splines'
#' for an interpolation with splines but has not been properly validated with simulations yet.
#' @param return.samples boolean, when set to TRUE the samples are returned instead of summary
#' statistics over the samples. Default is FALSE.
#' @param survival boolean, when set to TRUE the summary statistics over survival functions are
#' computed in addition to the summary statistics over the risk functions.
#' @param CIF boolean, when set to TRUE the summary statistics over cumulative incidence functions are
#' computed in addition to the summary statistics over the risk functions. Only applies to competing risks.
#' @param inv.link boolean, when set to TRUE the summary statistics are computed over the predictions of
#' longitudinal components after applying the inverse link function for each samples in addition to the
#' summary statistics over the linear predictors.
#' @param idLoop boolean, when set to TRUE the predictions are computed for each individual separately in a
#' loop, default is FALSE (all individuals in one call). This could be useful in case of very big prediction
#' tasks over low RAM devices, to avoid reaching limits.
#' @param resErrLong boolean, when set to TRUE the residual error for Gaussian or lognormal longitudinal
#' outcomes is added to the uncertainty of predictions (default is FALSE which predicts the true underlying
#' value of the longitudinal marker, i.e., error-free).
#' #' @param silentMode a boolean that will stop printing messages during computations if turned to TRUE.
#' @param ... Extra arguments.
#' @export
#' @importFrom Matrix bdiag Diagonal
#' @importFrom methods new

predict.INLAjoint <- function(object, newData=NULL, newDataSurv=NULL, timePoints=NULL, NtimePoints=50,
                              NsampleHY=20, NsampleFE=20, NsampleRE=50, id=NULL, Csurv=NULL, startTime=NULL,
                              horizon=NULL, baselineHaz="interpolation", return.samples=FALSE, FEonly=FALSE,
                              survival=FALSE, CIF=FALSE, inv.link=FALSE, idLoop=FALSE, resErrLong=FALSE, silentMode=FALSE, ...){
  # idGroup: loop over groups over random effects (useful if scaling issues)
  arguments <- list(...)
  # id is the id column name in dataset for survival data only (otherwise it's given by longitudinal)
  # Csurv is to get predictions conditional on survival up to given time
  ## strategy: dense (full sampling) ; update ; analytical
  # strategy =
  # joint: sample fixed effects and hyperparameters and for each sample,
  # estimate individual random effects posteriors and baseline risk
  # full: samples fixed effects and hyperparameters and for each sample, estimate
  # individual random effects posteriors (all in 1 call.
  # split: quantifies uncertainty from random effects conditional on the mode of hyperparameters
  # and adds it to the uncertainty from fixed effecst and hyperparameters.
  # (more stable and scalable compared to dense for same results (verified in sim))
  # dense: loops over samples of fixed effects and hyperparameters; for each sample,
  # estimate random effects individual posteriors for all new individuals
  # with one unique inla() call.
  # generic: loops over individuals; estimate all the random effects individual posteriors
  # for an individual for samples of fixed effects and hyperparameters with an inla generic0 call
  # this strategy does not account for survival time to compute the random effects posteriors
  # it is therefore more appropriate for longitudinal-only models.
  strategy="full3"
  Nsample=30
  NsampleREint=NULL
  SMPloop=FALSE
  if(is.null(NsampleREint) | strategy=="full" | strategy=="full2" | strategy=="full3" | strategy=="joint") NsampleREint <- Nsample
  if(is.null(newData)){ # if no new data is provided, return predicted fitted values
    PRED <- object$summary.fitted.values
    OUtc <- as.data.frame(object$.args$data$Yjoint)
    PRED$Outcome <- sapply(1:dim(PRED)[1], function(x) colnames(OUtc)[which(!is.na(OUtc[x,]))])
    return(PRED)
  }
  loopRE <- FALSE
  # boolean with default to FALSE. When 'NsampleRE' and 'Nsample' are large, the amount of
  # information to store in the random access memory of the computer can be large (creation of large matrices
  # for the computation of predictions), turning this boolean to TRUE will decompose the computation of
  # predictions to avoid reaching the limit of RAM of the computer (which would crash the program).
  if (!"INLAjoint" %in% class(object)){
    stop("Please provide an object of class 'INLAjoint' (obtained with joint() function).\n")
  }
  # baselineHaz = "smooth" | "interpolation"
  out <- NULL
  SumStats <- function(x) return(c(mean(x), sd(x), quantile(x, c(0.025, 0.5, 0.975))))
  if(!is.null(object$id)) id <- object$id else if(is.null(id)) stop("Please specify individual id column name with argument 'id'")
  is_Long <- is_Surv <- FALSE
  if(TRUE %in% c(unique(as.integer(newData[,object$id])) != 1:length(unique(as.integer(newData[,object$id]))))){
    # stop("Please use id = 1 to N individuals.")
    newID <- cbind(unique(as.integer(newData[,object$id])), 1:length(unique(as.integer(newData[,object$id]))))
    newData[,object$id] <- newID[sapply(newData[,object$id], function(x) which(newID[,1]==x)), 2]
    ChangeID <- TRUE
  }else{
    ChangeID <- FALSE
  }

  idVect <- na.omit(unique(object$.args$data[[paste0("ID", object[["REstruc"]][[1]])]]))
  if(!as.character(object["REstruc"])=="NULL"){
    is_Long <- TRUE
  }
  if(!is.null(object$SurvInfo)){
    if(is.null(idVect)){
      idVect <- unique(object$.args$data$expand1..coxph)
    }else{
      if(!any(idVect %in% unique(object$.args$data$expand1..coxph))) stop("id mismatch between longi and surv.")
    }
    is_Surv <- TRUE
    M <- length(object$survOutcome) # number of survival outcomes
    # check if we have baseline info for horizon
    for(m in 1:M){
      if(object$basRisk[[m]] %in% c("rw1", "rw2")){
        # add " method is interpolatiioon and forecast is required => switching to smooth
        # explain that interpolation means constant after last time point and suggest to use smooth to use a smooth prediction of baseline after
        if(horizon>max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]) & baselineHaz=="interpolation"){
          warning(paste0("The fitted model has baseline risk information up until value ",
                         max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]), " for survival outcome ", m, ". Since you ask for prediction at horizon ", horizon, " I will assume constant baseline hazard beyond the maximum available value. Alternatively, you can use baselineHaz='smooth' to use splines to predict the baseline hazard (for each sample). Alternatively, adding 'horizon' in the control options of the inla() call allows to extend the baseline beyond the last observed event time (linear extension based on last 2 values)."))
        }
      }
    }
  }
  # if(!is_Surv & strategy %in% c("dense", "split")){
  #   warning("strategy='generic' may be more efficient for this model.")
  #   # strategy="generic"
  # }
  if(!is_Long & is_Surv){
    if(exists("newDataSurv") & is.null(newData)){
      newData <- newDataSurv
    }else if(exists("newData") & is.null(newDataSurv)){
      newDataSurv <- newData
    }
  }
  if (inherits(newData, "tbl_df") || inherits(newData, "tbl")) {
    newData <- as.data.frame(newData)
  }
  if(!is_Long & !is_Surv) stop("Error, cannot recover ids from fitted model...")
  if(is_Surv & is.null(horizon)) stop("Please provide time horizon for prediction.")
  predL <- NULL
  predS <- NULL
  newPredS <- NULL
  if(is.null(object$id) & !is.null(id)) object$id <- id
  if(is.null(object$id)) stop("Please provide 'id' argument for new data.")
  ct <- object$misc$configs$contents
  if(is.null(ct)) stop("Please add argument 'cfg=TRUE' in control options when fitting the INLAjoint model to enable predictions.")
  if (ct$tag[1] == "Predictor") {
    ct$tag <- ct$tag[-1]
    ct$start <- ct$start[-1] - ct$start[2] + 1
    ct$length <- ct$length[-1]
  }
  paramVal <- object$misc$configs$config[[1]]$improved.mean # parameters value
  if(!is.null(startTime)) sTime <- startTime else sTime <- 0
  if(is.null(timePoints)){
    if(is.null(horizon)){
      timePoints <- seq(sTime, max(newData[, object$timeVar]), len=NtimePoints)
    }else{#} if(Csurv==0){
      timePoints <- seq(sTime, horizon, len=NtimePoints)
      # }else{
      #need to have a time point at Csurv there
    }
  }
  firstID <- unique(newData[, object$id])[1]
  if(strategy %in% c("generic", "dense", "split", "full", "full2", "full3", "joint")){
    if(!silentMode) message("Sample...")
    # set.seed(i)
    # SMPHt <- INLA::inla.hyperpar.sample(1, object)
    # SMPt <- INLA::inla.rjmarginal(1, object)
    # SMPH <- SMPHt
    # SMP <- SMPt
    # for(i in 1:(Nsample-1)){
    #   SMPH <- rbind(SMPH, SMPH[1])
    #   SMP[[1]] <- cbind(SMP[[1]], SMP[[1]][,1])
    # }
    if(strategy=="full2"){
      SMPH <- INLA::inla.hyperpar.sample(Nsample, object)
      SMPF <- INLA::inla.rjmarginal(Nsample*NsampleRE, object)
      SMP <- SMPF
      SMP$samples <- SMP$samples[, 1:Nsample]
    }else if(strategy=="full3"){
      SMPH <- INLA::inla.hyperpar.sample(NsampleHY, object)[rep(1:NsampleHY, each=NsampleFE),]
      SMP <- INLA::inla.rjmarginal(NsampleHY*NsampleFE, object)
      Nsample <- NsampleREint <- NsampleHY*NsampleFE
    }else{ # sample baseline for each random effects realizations
      SMPH <- INLA::inla.hyperpar.sample(Nsample, object)
      SMP <- INLA::inla.rjmarginal(Nsample, object)
    }
    # if(Nsample==10){
    #   load("~/Documents/GitHub/INLAjoint/SMP10.RData")
    #
    # }else if(Nsample==11){
    #   load("~/Documents/GitHub/INLAjoint/SMP11.RData")
    #
    # }
    if(FEonly){
      for(m in 1:M){
        idB_H <- c(sapply(1:M, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                            (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                               ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))
        SMP$samples[idB_H,] <- rowMeans(SMP$samples[idB_H,])
      }
    }
    if(!silentMode) message(paste0("Compute predictions..."))
  }
  ParVal <- new("dgTMatrix", Dim=c(sum(ct$length), as.integer(Nsample)))
  if(strategy %in% c("dense", "split", "full", "full2", "full3", "joint")){ # use mode of hyperpar to estimate random effects
    ParValMode <- object$misc$configs$config[[1]]$improved.mean
  }
  if(is_Surv){
    RWBH <- which(object$basRisk %in% c("rw1", "rw2"))
    if(length(RWBH)>0){
      SMP$samples[,1] <- ParValMode[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                                 (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                                    ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
                                      ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                       !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      if(strategy=="full2"){
        SMPF$samples[,1] <- ParValMode[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                                    (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                                       ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
                                         ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                          !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      }
      ParVal[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                          (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                             ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
               ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                !ct$tag %in% paste0("baseline", 1:M, ".hazard"))]),] <- SMP$samples#[which(!substr(rownames(SMP$samples),1, 9) %in% paste0("baseline", 1:M)),]
      # remove samples for now to check estimation if it works (REMOVE THIS AFTER - USELESS)
      # ParVal[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
      #                     (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
      #                        ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- NA

    }else{
      SMP$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                       !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      if(strategy=="full2"){
        SMPF$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                         !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      }
      ParVal[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                !ct$tag %in% paste0("baseline", 1:M, ".hazard"))]),] <- SMP$samples#[which(!substr(rownames(SMP$samples),1, 9) %in% paste0("baseline", 1:M)),]
    }
    # if(baselineHaz=="smooth"){
    #   for(m in 1:M){
    #     BHpos <-which(ct$tag %in% paste0("baseline", m, ".hazard"))
    #     smoothBL <- function(x){
    #       splBH <- splinefun(object$summary.random[[paste0("baseline", m, ".hazard")]]$ID,
    #                          x[c(ct$start[BHpos]:(ct$start[BHpos]+ct$length[BHpos]-1))])
    #       splBH(TPO)
    #     }
    #     NewBas <- apply(SMP$samples, 2, smoothBL)
    #     if(ct$start[BHpos]==1){
    #       ParVal <- INLA::inla.as.sparse(rbind(NewBas, ParVal[-c(ct$start[BHpos]:(ct$start[BHpos]+ct$length[BHpos]-1)),]))
    #     }else if(ct$start[BHpos]>1){
    #       ParVal <- INLA::inla.as.sparse(rbind(ParVal[1:(ct$start[BHpos]-1),],
    #                                      NewBas, ParVal[-c(1:(ct$start[BHpos]+ct$length[BHpos]-1)),]))
    #     }
    #     ct$start[-c(1:BHpos)] <- ct$start[-c(1:BHpos)] - ct$length[BHpos] + length(TPO)
    #     ct$length[BHpos] <- length(TPO)
    #   }
    # }
  }else{
    # use mode as first sample
    SMP$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))])]
    if(strategy=="full2"){
      SMPF$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))])]
    }
    ParVal[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))]),] <- SMP$samples
  }
  nRE <- 0
  if(is_Long | !is.null(object[["REstrucS"]])){
    K <- length(object$famLongi) # number of longitudinal outcomes
    ct2 <- ct
    lenPV <- length(paramVal)
    SMPsel <- which(ct$length==1 &
                      substr(ct$tag, nchar(ct$tag)-2, nchar(ct$tag)-1)=="_L" |
                      substr(ct$tag, nchar(ct$tag)-3, nchar(ct$tag)-2)=="_L") # if >10 markers
    NamesH <- colnames(SMPH)
    nRE <- length(object[["REstruc"]])
    if(is.null(object[["REstrucS"]])){
      if(nRE==1){
        BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample))) # adapt size
        diag(BD_Cmat) <- sqrt(1/SMPH[, which(substr(colnames(SMPH), 1, 16)=="Precision for ID")])
      }else if(nRE>1){
        # identify the position of the cholesky elements in hyperparameters
        if(object$corLong){










































          PosH <- which(substr(NamesH, 1, 5)=="Theta" &
                          substr(NamesH, nchar(NamesH)-nchar(object[["REstruc"]][1])-1,
                                 nchar(NamesH))==paste0("ID", object[["REstruc"]][1]))
          if(strategy%in%c("generic", "dense", "split", "full", "full2", "full3", "joint")){
            # Block-Diagonal Cmatrix for all samples
            BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample))) # adapt size
            # make function that compute the precision matrix and place it in BD_Cmat



            L <- matrix(0, nrow=nRE, ncol=nRE)
            # function to convert cholesky to precision
            Chol_Prec <- function(x){
              diag(L) <- exp(x[PosH][1:nRE])
              L[lower.tri(L)] <- x[PosH][-c(1:nRE)]
              return(L %*% t(L))
            }
            SMP_prec <- apply(SMPH, 1, Chol_Prec)
            # indices for BC_Cmat
            ind_BD_Cmat <- cbind(rep(1:(Nsample*nRE), each=nRE), rep(1:nRE, Nsample)+rep(nRE*(1:Nsample-1), each=nRE^2))
            # fill BD_Cmat
            BD_Cmat[ind_BD_Cmat] <- c(SMP_prec)
          }else if(strategy==20){
            #sobol_init = sobol(n=1, dim=nRE, init = TRUE, scrambling = 1, normal = TRUE)
            # sobol_sample = sobol(NsampleRE, dim = nRE, scrambling = 1, init = FALSE, normal = TRUE)
            # chol_re <- matrix(0, nrow = nRE, ncol=nRE)
            # diag(chol_re) <- exp(SMPH[1,PosH][1:nRE])
            # chol_re[lower.tri(chol_re)] <- SMPH[1,PosH][-c(1:nRE)]
            #
            # str(SMPH)
            # solve(chol_re) %*% sobol_sample[1,]
            #
            # A <- SMPH[1,PosH]
            # B <- sobol_sample[1,]

            # RE_values2 <- mvtnorm::rmvnorm(NsampleRE, sigma=solve(Cmatrix))
            # RE_values <-t(sapply(1:nRE, function(x) RE_values2[, seq(x, Nsample*nRE, by=nRE)]))
            # need to weight samples with probability density from observations!
          }
        }else{
          nRE_pk <- 1
          # Block-Diagonal Cmatrix for all samples
          BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample), as.integer(nRE*Nsample))) # adapt size
          for(k in 1:K){
            PosH <- which(substr(NamesH, 1, 5)=="Theta" &
                            substr(NamesH, nchar(NamesH)-nchar(object[["REstruc"]][nRE_pk])-1,
                                   nchar(NamesH))==paste0("ID", object[["REstruc"]][nRE_pk]))
            nRE_k <- length(which(substr(object[["REstruc"]], nchar(object[["REstruc"]])-2, nchar(object[["REstruc"]]))==paste0("_L", k) |
                                    substr(object[["REstruc"]], nchar(object[["REstruc"]])-3, nchar(object[["REstruc"]]))==paste0("_L", k)))
            if(nRE_k==1){
              SMP_prec_k <- sqrt(1/SMPH[,which(substr(colnames(SMPH), 1, 16)=="Precision for ID" &
                                                 (substr(colnames(SMPH), nchar(colnames(SMPH))-2, nchar(colnames(SMPH)))==paste0("_L", k) |
                                                    substr(colnames(SMPH), nchar(colnames(SMPH))-3, nchar(colnames(SMPH)))==paste0("_L", k)))])
            }else{
              L <- matrix(0, nrow=nRE_k, ncol=nRE_k)
              # function to convert cholesky to precision
              Chol_Prec <- function(x){
                diag(L) <- exp(x[PosH][1:nRE_k])
                L[lower.tri(L)] <- x[PosH][-c(1:nRE_k)]
                return(L %*% t(L))
              }
              SMP_prec_k <- apply(SMPH, 1, Chol_Prec)
            }
            # indices for BC_Cmat
            ind_BD_Cmat_k <- cbind(rep(rep(1:nRE_k, each=nRE_k), Nsample)+(rep(seq(nRE_pk, (Nsample*nRE), by=nRE), each=nRE_k^2)-1),
                                   rep(1:nRE_k, Nsample*nRE_k)+(rep(seq(nRE_pk, (Nsample*nRE), by=nRE), each=nRE_k^2)-1))
            # fill BD_Cmat
            BD_Cmat[ind_BD_Cmat_k] <- c(SMP_prec_k)
            nRE_pk <- nRE_pk + nRE_k # go to next block
          }
        }
      }
    }else{ # if there is at least a frailty, need to do the full model
      if(is_Long) nRES <- length(object[["REstrucS"]]) else nRES <- nRE
      if(nRE==1){ # only frailty
        BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample))) # adapt size
        diag(BD_Cmat) <- sqrt(1/SMPH[, which(substr(colnames(SMPH), 1, 16)=="Precision for ID")])
      }else if(nRE>1){
        # need to do the full model to get posteriors to sample from as there is at least longi and frailty here

      }



















    }
    # if(strategy %in% c("dense", "split", "full", "joint")){
    #   RE_values2 <- mvtnorm::rmvnorm(NsampleRE, sigma=solve(BD_Cmat))
    #   RE_values <-t(sapply(1:nRE, function(x) RE_values2[, seq(x, Nsample*nRE, by=nRE)]))
    #   # need to weight samples with probability density from observations!
    # }
    # RE_values2 <- mvtnorm::rmvnorm(NsampleRE, sigma=solve(BD_Cmat))
    # RE_values <-t(sapply(1:nRE, function(x) RE_values2[, seq(x, Nsample*nRE, by=nRE)]))
    ResErrFixed <- vector("list", K)
    if(is.null(object[["REstrucS"]])){
      REnames <- c(sapply(object["REstruc"], function(x) paste0("ID", x)))
    }else{
      if(!as.character(object["REstruc"])=="NULL"){
        REnames <- c(sapply(object["REstruc"], function(x) paste0("ID", x)))
        REnamesS <- object[["REstrucS"]]
      }else{
        REnames <- REnamesS <- object[["REstrucS"]]
      }
    }
    posRE <- ct$start[sapply(REnames, function(x) which(ct$tag==x))]
    assocNs <- object$assoc
    assocNa <- object$assoc_Names
    if(is_Surv){
      assocPos <- sapply(assocNs, function(x) grep(x, ct$tag))
      # identify the longitudinal needed for association
      # first identify shared part from longitudinal (no duplicates, so if CV from longitudinal 1 is shared twice, we need to repeat it)
      # outcomeAssoc <- names(object$.args$data$Yjoint)[which(substr(names(object$.args$data$Yjoint), 1, nchar(names(object$.args$data$Yjoint))-1) %in% ct2$tag)]


      OutcNam <- substr(names(object$.args$data$Yjoint), 1, nchar(names(object$.args$data$Yjoint))-1)
      outcomeAssoc <- names(object$.args$data$Yjoint)[unlist(sapply(1:length(OutcNam), function(x) if(length(grep(OutcNam[x], ct2$tag))!=0) return(x)))]
      outcomeAssoc2 <- sapply(outcomeAssoc, function(x) strsplit(x, split = "_S")[[1]][1])
      requiredAssoc <- sapply(assocNs, function(x) strsplit(x, split = "_S")[[1]][1])
      patternAsso <- unname(sapply(requiredAssoc, function(x) which(x==outcomeAssoc2)))
      # REnames <- object[["REstruc"]]#c(substr(object[["REstrucS"]], 3, nchar(object[["REstrucS"]])),
      if(length(object[["REstrucS"]])>0){
        FRAIL_ind <- 1:length(object[["REstrucS"]])
      }else{
        FRAIL_ind <- NULL
      }
      SRE_indi <- NULL
      if(!is.null(REnames)){
        for(i in 1:length(REnames)){
          if(length(grep(REnames[i], assocNs))>0){
            SRE_indi <- c(SRE_indi, length(object[["REstrucS"]])+i)
          }
        }
        SRE_inda <- unlist(sapply(REnames, function(x) grep(x, assocNs)))# indicator for associations that are SRE_ind
      }else{
        SRE_inda <- NULL
      }
      # SRE_indi <- which(sapply(patternAsso, function(x) length(x))==0) # identify SRE_ind
      if(length(SRE_inda)>0){
        # SRE_indipos <- assocPos[SRE_inda]
        # SRE_indiN <- assocNs[SRE_inda]
        patternAsso2 <- unname(unlist(patternAsso[-SRE_inda]))
        patternAsso <- 1:length(patternAsso)
        assocNs2 <- assocNs[-SRE_inda]
      }else{
        patternAsso2 <- patternAsso
        assocNs2 <- assocNs
      }
    }
  }
  # if(length(unique(newData[, object$id]))==1) Nthreads=1
  # if(is.null(Nthreads)) Nthreads = future::availableCores("mc.cores")-1
  # registerDoFuture()
  # registerDoRNG()
  # message(paste0("Number of threads: ", Nthreads))
  # plan("multisession", workers = Nthreads)
  # PRED <- foreach(idPred = unique(newData[, object$id]), .export=c(ls(), ls(envir=.GlobalEnv), lsf.str(envir=.GlobalEnv), "joint")) %dopar% {
  if(idLoop) idLoopSet <- FALSE else idLoopSet <- TRUE # used when all individuals random effects done in 1 call
  for(idPred in unique(newData[, object$id])){
    if(idLoop | strategy=="generic"){
      ND <- newData[newData[, object$id] == idPred,] # select only one id at a time
    }else{
      #IDP_ = unique(newData[, object$id])[length(unique(newData[, object$id]))] # skip the loop
      ND <- newData
    }
    if(!is.null(object$lonFacChar) & length(which(names(object$lonFacChar) %in% colnames(ND)))>0){
      for(Fi in which(names(object$lonFacChar) %in% colnames(ND))){
        # colClass <- apply(ND, 2, class)
        ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])] <- factor(gsub("[^[:alnum:] ]","", ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])]), levels=gsub("[^[:alnum:] ]","", object$lonFacChar[[Fi]]))
        ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])] <- factor(gsub(" ","", ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])]), levels=gsub(" ","", object$lonFacChar[[Fi]]))
        # ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])] <- factor(ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])], levels=object$lonFacChar[[Fi]])
      }
    }
    if(is_Long & strategy %in% c("dense", "split", "full", "full2", "full3", "joint")){ # add time points at observed longi to compute density
      if(is.null(Csurv)){
        TPO <- sort(unique(c(timePoints, ND[, object$timeVar])))
        NTP <- length(TPO)
      }else if(Csurv==0){
        TPO <- sort(unique(c(timePoints, ND[, object$timeVar])))
        NTP <- length(TPO)
      }else{
        TPO <- sort(unique(c(timePoints, Csurv, ND[, object$timeVar])))
        NTP <- length(TPO)
      }
    }else if(is_Long & is.null(Csurv)){
      TPO <- sort(unique(c(timePoints, max(ND[, object$timeVar]))))
      NTP <- length(TPO)
    }else if(!is_Long & is.null(Csurv)){
      TPO <- timePoints
      NTP <- NtimePoints
    }else if(Csurv==0){
      TPO <- timePoints
      NTP <- NtimePoints
    }else{
      TPO <- sort(c(timePoints, Csurv))
      NTP <- NtimePoints+1
    }
    call.new2 <- object$call
    TXT1 <- NULL
    if(is_Surv){
      if(is_Long & !is.null(newDataSurv)){
        NDS <- newDataSurv
      }else{
        if(!is.null(object[["REstrucS"]])){
          NDS <- ND
        }else{
          NDS <- ND[which(!duplicated(ND[,id], fromLast=TRUE)),]
        }
        for(si in 1:length(object$SurvInfo)){
          if(length(as.character(object$SurvInfo[[si]]$survOutcome))==1){
            S_Outc <- as.character(object$SurvInfo[[si]]$survOutcome)
          }else if(length(as.character(object$SurvInfo[[si]]$survOutcome))==3){
            S_Outc <- as.character(object$SurvInfo[[si]]$survOutcome)[3]
          }
          if(length(as.character(object$SurvInfo[[si]]$nameTimeSurv))==1){
            S_nam <- as.character(object$SurvInfo[[si]]$nameTimeSurv)
            T_nam <- as.character(object$SurvInfo[[si]]$nameTrunc)
          }else if(length(as.character(object$SurvInfo[[si]]$nameTimeSurv))==3){
            S_nam <- as.character(object$SurvInfo[[si]]$nameTimeSurv)[3]
            T_nam <- as.character(object$SurvInfo[[si]]$nameTrunc)[3]
          }
          if(!S_Outc %in% colnames(NDS)){
            NDS <- cbind(NDS, 0)
            colnames(NDS)[length(colnames(NDS))] <- S_Outc
          }
          if(!S_nam %in% colnames(NDS)){
            if(is.null(object$timeVar)){
              colTS <- which(colnames(ND) %in% unlist(sapply(object$SurvInfo, function(x) x$nameTimeSurv)))
              mTS <- max(ND[colTS])
            }else{
              mTS <- max(ND[object$timeVar])
            }
            NDS <- cbind(NDS, mTS)
            colnames(NDS)[length(colnames(NDS))] <- S_nam
          }else{
            # set time > 0 if only providing longitudinal measurements at time 0?
          }
          if(length(T_nam)>0){
            if(!is.na(T_nam)){
              if(!T_nam %in% colnames(NDS)){
                NDS <- cbind(NDS, 0)
                colnames(NDS)[length(colnames(NDS))] <- T_nam
              }
            }
          }
        }
      }
      if(max(ND[object$timeVar])>horizon & is.null(Csurv)) warning(paste0("horizon = ", horizon, " and there are observations up to ", max(ND[object$timeVar]), ". It is likely not what you want but you can use Csurv argument if you want to to force predictions conditional on future observations."))
      horizon2 <- max(TPO)+0.0001
      # SdataPred <- ND[!duplicated(ND[, object$id]),]
      if(!is.null(object[["REstrucS"]])){
        SdataPred <- NDS[NDS[,id]==idPred,]
      }else{
        SdataPred <- NDS[NDS[,id]==idPred,][1,]
      }
      #if(!is.null(object$dataSurv))
      for(m in 1:M){
        if(length(paste0(object$SurvInfo[[m]]$nameTimeSurv))>1) NTS <- tail(paste0(object$SurvInfo[[m]]$nameTimeSurv),1) else NTS <- paste0(object$SurvInfo[[m]]$nameTimeSurv)
        if(length(paste0(object$SurvInfo[[m]]$survOutcome))>1) SVO <- tail(paste0(object$SurvInfo[[m]]$survOutcome),1) else SVO <- paste0(object$SurvInfo[[m]]$survOutcome)
        if(NTS %in% colnames(SdataPred)){
          SdataPred[nrow(SdataPred), NTS] <- horizon2
        }else{
          SdataPred <- cbind(SdataPred, horizon2)
          colnames(SdataPred)[length(colnames(SdataPred))] <- NTS
        }
        if(!(SVO %in% colnames(SdataPred))){
          SdataPred$SVO <- 0
          colnames(SdataPred)[which(colnames(SdataPred)=="SVO")] <- SVO
        }else{
          if(SdataPred[nrow(SdataPred), SVO]==1 & !is.null(object[["REstrucS"]])){
            SdataPred <- SdataPred[c(1:nrow(SdataPred), nrow(SdataPred)),]
          }
          SdataPred[nrow(SdataPred), SVO] <- 0
        }
        # remove truncation to force predict at all given time points
        if(!is.null(object$SurvInfo[[m]]$nameTrunc) & !is.null(startTime)){
          SdataPred[, which(colnames(SdataPred)==object$SurvInfo[[m]]$nameTrunc)] <- min(TPO)
        }
      }
      if(!is.null(object$dataSurv)){
        if(paste0(object$dataSurv)[1]=="list"){
          if(length(object[["REstrucS"]])>9) stop("Predictions not implemented for 10+ frailties, contact INLAjoint@gmail.com")
          for(m in 1:(length(paste0(object$dataSurv))-1)){ # all lines
            if(length(grep(paste0("_S", m), substr(object[["REstrucS"]],
                                            start=nchar(object[["REstrucS"]])-2,
                                            stop=nchar(object[["REstrucS"]]))))>0){
              assign(paste0(object$dataSurv)[m+1], SdataPred)
            }else{ # only last line
              assign(paste0(object$dataSurv)[m+1], SdataPred[nrow(SdataPred),])
            }
          }
        }else{
          if(length(grep("_S1", substr(object[["REstrucS"]],
                                                 start=nchar(object[["REstrucS"]])-2,
                                                 stop=nchar(object[["REstrucS"]]))))>0){
            assign(paste0(object$dataSurv)[m+1], SdataPred)
          }else{ # only last line
            assign(paste0(object$dataSurv), SdataPred[nrow(SdataPred),])
          }
          assign(paste0(object$dataSurv), SdataPred)
        }
      }else{
        object$dataSurv <- SdataPred
      }
      CTP <- c(TPO, max(TPO)+tail(diff(TPO),1)) # add fake point to extend and have the last value
      # if(!is.null(object$cutpoints)) assign(paste0(object$cutpoints), CTP) else TXT1 <- ", cutpoints=CTP"
    }else{
      NDS <- NULL
      M <- 0
    }

    if(is_Long){
      LdataPred <- ND[rep(1, length(TPO)), ]
      LdataPred[, object$timeVar] <- TPO
      # if(!is.null(object$lonFacChar) & length(which(names(object$lonFacChar) %in% colnames(LdataPred)))>0){
      #   for(Fi in which(names(object$lonFacChar) %in% colnames(LdataPred))){
      #     # colClass <- apply(LdataPred, 2, class)
      #     LdataPred[, which(colnames(LdataPred)==names(object$lonFacChar)[Fi])] <- factor(LdataPred[, which(colnames(LdataPred)==names(object$lonFacChar)[Fi])], levels=object$lonFacChar[[Fi]])
      #   }
      # }
      # if(!is.null(object$survFacChar) & length(which(names(object$survFacChar) %in% colnames(SdataPred)))>0){
      #   for(Fi in which(names(object$survFacChar) %in% colnames(SdataPred))){
      #     # colClass <- apply(SdataPred, 2, class)
      #     SdataPred[, which(colnames(SdataPred)==names(object$survFacChar)[Fi])] <- factor(SdataPred[, which(colnames(SdataPred)==names(object$survFacChar)[Fi])], levels=object$survFacChar[[Fi]])
      #   }
      # }
      if(is.null(startTime)){ # start predictions from first observed time
        TPO <- TPO[which(TPO>=min(ND[which(ND[, object$id]==idPred), object$timeVar]))]
        NTP <- length(TPO)
        LdataPred <- LdataPred[which(LdataPred[, object$timeVar]>=min(ND[which(ND[, object$id]==idPred), object$timeVar])),]
      }
      if(!is.null(object$dataLong)){
        if(paste0(object$dataLong)[1]=="list"){
          for(m in 1:(length(paste0(object$dataLong))-1)){
            assign(paste0(object$dataLong)[m+1], LdataPred)
          }
        }else{
          assign(paste0(object$dataLong), LdataPred)
        }
      }else{
        object$dataLong = LdataPred
      }
      # if(!is.null(object$dataLong)) assign(paste0(object$dataLong), LdataPred)
      # if(!is.null(object$dataSurv) & is_Surv) assign(paste0(object$dataSurv), SdataPred)
    }
    if(length(grep("dataSurv", object$call))==0 & is_Surv){
      call.new2[[length(object$call)]] <- paste(c(substr(object$call[[length(object$call)]],
                                                         start=1,
                                                         stop=nchar(object$call[[length(object$call)]])-1),
                                                  ", dataSurv=SdataPred, dataOnly=TRUE)"), collapse='')
    }else{
      call.new2[[length(object$call)]] <- paste(c(substr(object$call[[length(object$call)]],
                                                         start=1,
                                                         stop=nchar(object$call[[length(object$call)]])-1),
                                                  ", dataOnly=TRUE)"), collapse='')
    }
    call.new2 <- paste(call.new2, collapse='')
    if(is_Surv){
      # fix warning when assigning ctp to a vector of numerics instead of a name
      if(length(grep("control = list\\(", call.new2))>0){
        call.new2 <- gsub("control = list\\(", "control = list\\(cutpointsF = CTP,", call.new2)
      }else{
        call.new2 <- paste0(substr(call.new2,
                                   start=1,
                                   stop=nchar(call.new2)-1),
                            ", control=list(cutpointsF=CTP))")
      }
    }
    NEWdata <- suppressWarnings(eval(parse(text=call.new2))) # maybe need to store functions of time in the object?
    survPart <- NULL
    # if(is_Surv) survPart <- unlist(sapply(1:M, function(x) which(!is.na(eval(parse(text=paste0("NEWdata$baseline", x, ".hazard.idx")))))))
    if(is_Surv) survPart <- c(unlist(sapply(1:M, function(x) which(!is.na(eval(parse(text=paste0("NEWdata$y", x, "..coxph"))))))))
    if(!is.list(NEWdata)) NEWdata <- as.list(as.data.frame(NEWdata))
    if(is_Long | !is.null(object[["REstrucS"]]) | strategy=="joint"){
      ###              ###
      ### LONGITUDINAL ###
      ###              ###
      ctL <- ct
      if(F){#idPred %in% idVect & identical(ND[, object$timeVar], object$.args$data[[paste0(object$timeVar, "_L1")]][which(object$.args$data$IDIntercept_L1==idPred)])){
        # no new data, we can use the current model as it is to predict over the timePoints
        LdataPred <- ND[rep(1, length(TPO)), ]
        LdataPred[, object$timeVar] <- TPO
        # if(is_Long) predL <- rbind(predL, predLongi(object, LdataPred))
      }else if (idPred %in% idVect & !(identical(ND[, object$timeVar], object$.args$data[[paste0(object$timeVar, "_L1")]][which(object$.args$data$IDIntercept_L1==idPred)])) & strategy==3){
        # if id already known (random effects available) and new observations and update method
        # new data available, we need to update the model
        # first update the data with new observations
        # if(is.null(oldData)) stop("Please provide original dataset 'oldData' for the 'update' strategy.")
        # ND2 <- rbind(oldData, ND)
        # ND2 <- ND2[!duplicated(ND2),]
        # ND2 <- ND2[order(ND2[, object$id], ND2[, object$timeVar]),]
        # assign(paste0(object$dataLong), ND2)
        # # create data with dataOnly=T and replace in object to rerun.
        # call.new <- object$call
        # call.new[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
        #                                               start=1,
        #                                               stop=nchar(object$call[[length(object$call)]])-2),
        #                                        substr(object$call[[length(object$call)]],
        #                                               start=nchar(object$call[[length(object$call)]]),
        #                                               stop=nchar(object$call[[length(object$call)]])),
        #                                        ", dataOnly=TRUE)")
        # updatedData <- eval(parse(text=call.new))
        # object.new <- object
        # object.new$.args$data <- updatedData
        # object.new$.args$offset <- rep(0, length(object.new$.args$data[[1]]))
        # object.new$.args$E <- updatedData$E..coxph
        # print("Updating model with new data.")
        # object.new <- joint.rerun(object.new) # this is the updated fitted model
        # print("Computing longitudinal predictions.")
        # # now we can do the prediction
        # LdataPred <- ND[rep(1, length(TPO)), ]
        # LdataPred[, object.new$timeVar] <- TPO
        #
        # if(!is.null(object[["REstruc"]])){
        #   IDshift <- 0
        #   LengthUniqueID <- length(unique(na.omit(object$.args$data[[REnames[1]]])))
        #   if(object$corLong){
        #     for(RE in REnames){
        #       ctL$start[which(ctL$tag==RE)] <- ctL$start[which(ctL$tag==RE)] + as.integer(unique(LdataPred[,object$id])) + IDshift*LengthUniqueID -1
        #       IDshift <- IDshift + 1
        #     }
        #   }else{
        #     for(k in 1:K){
        #       IDshift <- 0
        #       Ndigits <- ifelse(k<10, 0, 1)
        #       RE_k <- which(substr(ctL$tag, 1, 2)=="ID" & substr(ctL$tag, nchar(ctL$tag)-Ndigits, nchar(ctL$tag))==k)
        #       for(RE in RE_k){
        #         ctL$start[RE] <- ctL$start[RE] + as.integer(unique(LdataPred[,object$id])) + IDshift*LengthUniqueID -1
        #         IDshift <- IDshift + 1
        #       }
        #     }
        #   }
        # }
        # if(!is.null(object$dataLong)) assign(paste0(object$dataLong), LdataPred)
        # call.new <- object$call
        # call.new[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
        #                                               start=1,
        #                                               stop=nchar(object$call[[length(object$call)]])-1),
        #                                        ", dataOnly=TRUE, longOnly=TRUE)")
        # NEWdata <- eval(parse(text=call.new))
        # NEWdata[REnames] <- NEWdata[paste0("W", object[["REstruc"]])]
        #
        # ANew <- matrix(0, nrow=length(NEWdata[[1]]), ncol=sum(ctL$length))
        # ANew[, ctL$start[which(ctL$tag %in% names(NEWdata))]] <- do.call(cbind, NEWdata[ctL$tag])
        # ANew <- INLA::inla.as.sparse(ANew, na.rm=TRUE)
        # # get sd from Q matrix
        # QQ= object.new$misc$configs$config[[1]]$Q
        # d.QQ=diag(as.matrix(QQ))
        # QQ=as.matrix(QQ) + t(as.matrix(QQ))
        # diag(QQ) = d.QQ
        # predL <- data.frame(rep(LdataPred[, object.new$id], K), rep(LdataPred[, object.new$timeVar], K),
        #                     rep(object.new$longOutcome, each=NTP), as.matrix(ANew %*% paramVal),
        #                     sqrt(diag(as.matrix(ANew) %*% solve(QQ, t(as.matrix(ANew))))))
        #
        # colnames(predL) <- c(object.new$id, object.new$timeVar, "Outcome", "Mean", "Sd")
        # return(predL)
      }else if (!(idPred %in% idVect) & !(identical(ND[, object$timeVar], object$.args$data[[paste0(object$timeVar, "_L1")]][which(object$.args$data$IDIntercept_L1==idPred)])) & strategy==3){
        stop("'update' strategy not available for prediction over new individuals (only for new observations of existing individuals), please select another strategy.")

      }else if(strategy %in% c("generic", "dense", "split", "full", "full2", "full3", "joint")){
        if(strategy=="full3"){
          inND <- as.integer(ND[, object$id])
          ND <- ND[rep(1:nrow(ND), NsampleFE),]
          ND[, object$id] <- rep(inND, NsampleFE) + rep(max(inND), NsampleFE*length(inND))*rep(0:(NsampleFE-1), each=length(inND))
          inNDS <- as.integer(NDS[, object$id])
          NDS <- NDS[rep(inNDS, NsampleFE),]
          NDS[, object$id] <- rep(inNDS, NsampleFE) + rep(max(inNDS), NsampleFE*length(inNDS))*rep(0:(NsampleFE-1), each=length(inNDS))
        }
        # convert data into INLAjoint format with dataOnly option
        if(!is.null(object$dataLong) & is_Long){
          if(paste0(object$dataLong)[1]=="list"){
            for(m in 1:(length(paste0(object$dataLong))-1)){
              assign(paste0(object$dataLong)[m+1], ND)
            }
          }else{
            assign(paste0(object$dataLong), ND)
          }
        }
        if(!is.null(object[["REstrucS"]]) | (strategy %in% c("dense", "split", "full", "full2", "full3", "joint") & is_Surv)){
          if(paste0(object$dataSurv)[1]=="list"){
            for(m in 1:(length(paste0(object$dataSurv))-1)){
              if(length(grep(paste0("_S", m), substr(object[["REstrucS"]],
                                                     start=nchar(object[["REstrucS"]])-2,
                                                     stop=nchar(object[["REstrucS"]]))))>0){
                assign(paste0(object$dataSurv)[m+1], NDS)
              }else{ # only last line
                assign(paste0(object$dataSurv)[m+1], NDS)
              }
            }
          }else{
            if(length(grep("_S1", substr(object[["REstrucS"]],
                                         start=nchar(object[["REstrucS"]])-2,
                                         stop=nchar(object[["REstrucS"]]))))>0){
              assign(paste0(object$dataSurv)[m+1], NDS)
            }else{ # only last line
              assign(paste0(object$dataSurv), NDS)
            }
          }
          call.new <- paste(object$call, collapse='')
          CTP <- object$.args$data$baseline1.hazard.values
          if(length(grep("control = list\\(", call.new)>0)){
            call.new <- gsub("control = list\\(", "control = list\\(cutpointsF = CTP,", call.new)
            call.new <- paste(substr(call.new,
                                     start=1,
                                     stop=nchar(call.new)-1),
                              ", dataOnly=TRUE)", collapse='')
          }else{
            call.new <- paste0(substr(call.new,
                                      start=1,
                                      stop=nchar(call.new)-1),
                               ", dataOnly=TRUE, control=list(cutpointsF=CTP))")
          }
          FRM <- object$.args$formula
          FRM2 <- paste0(paste0(FRM)[2], paste0(FRM[1], paste0(FRM[3])))
          SPLIT_n <- strsplit(FRM2, ", n = (.*?),")[[1]] # change the length of iid random effects as data is different
          # recover length of each iid random effect groups
          nre_pr <- NULL
          if(is_Long & is.null(object[["REstrucS"]])){
            if(object$corLong) rmvCL = 0 else rmvCL=1
            if(length(object$famLongi) != (length(SPLIT_n)-rmvCL)) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
            for(nre_p in 1:length(object$famLongi)){
              if(nre_p<10){
                nre_10p = 0
              }else if(nre_p>=10){
                nre_10p = 1
              }
              nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
            }
            if(object$corLong) nre_pr <- sum(nre_pr)
          }else if(is_Long & !is.null(object[["REstrucS"]])){
            if(object$corLong){
              if((1+length(object[["REstrucS"]])) != (length(SPLIT_n)-1)) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
            }else{
              if((length(object$famLongi)+length(object[["REstrucS"]])) != (length(SPLIT_n)-1)) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
            }
            for(nre_p in 1:length(object$famLongi)){ # first longi random effects
              if(nre_p<10){
                nre_10p = 0
              }else if(nre_p>=10){
                nre_10p = 1
              }
              nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]])))))
            }
            if(object$corLong) nre_pr <- sum(nre_pr)
            for(nre_p in 1:length(object[["REstrucS"]])){ # then surv frailty random effects
              if(nre_p<10){
                nre_10p = 0
              }else if(nre_p>=10){
                nre_10p = 1
              }
              nre_pr <- c(length(grep(paste0("_S", nre_p), substr(object[["REstrucS"]], start=nchar(object[["REstrucS"]])-2-nre_10p, stop=nchar(object[["REstrucS"]])))), nre_pr)
            }
          }else if(!is_Long & !is.null(object[["REstrucS"]])){
            if(length(object[["REstrucS"]]) != (length(SPLIT_n)-1)) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
            for(nre_p in 1:length(object[["REstrucS"]])){
              if(nre_p<10){
                nre_10p = 0
              }else if(nre_p>=10){
                nre_10p = 1
              }
              nre_pr <- c(nre_pr, length(grep(paste0("_S", nre_p), substr(object[["REstrucS"]], start=nchar(object[["REstrucS"]])-2-nre_10p, stop=nchar(object[["REstrucS"]])))))
            }
          }
          if(object$corLong){
            FRM3 <- paste(paste(sapply(1:(1+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], ", n = ", nre_pr[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
          }else{
            if((length(object$famLongi)+length(object[["REstrucS"]]))>0){
              FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], ", n = ", nre_pr[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
            }else{
              FRM3 <- FRM2
            }
          }
          # remove baseline from formula as it comes from full model (sampled as fixed effect)
          if(is_Surv & strategy!="joint"){
            nRWBH <- which(unlist(object$basRisk) %in% c("rw1", "rw2"))
            if(length(nRWBH)>0){
              for(mBH in nRWBH){
                FRM4 <- FRM3
                SPLIT_BH <- strsplit(FRM4, paste0("[:+:] f\\(baseline", mBH, ".hazard"))[[1]][1]
                SPLIT_BH2 <- strsplit(FRM4, "constr = FALSE, diagonal = 0.01, scale.model = TRUE\\)")[[1]][-1]
                if(length(SPLIT_BH2)==1){
                  FRM3 <- paste0(SPLIT_BH, SPLIT_BH2)
                }else if(length(SPLIT_BH2)>1){
                  for(BHi in 1:(length(SPLIT_BH2)-1)){
                    SPLIT_BH3 <- SPLIT_BH2
                    SPLIT_BH2 <- paste(paste0(SPLIT_BH3[1], "constr = FALSE, diagonal = 0.01, scale.model = TRUE)"), SPLIT_BH3[-1], collapse='')
                  }
                }
                FRM3 <- paste(SPLIT_BH, SPLIT_BH2, collapse='')
              }
            }
          }
        }else{
          FRM <- object$.args$formula
          FRM2 <- paste0(paste0(FRM)[2], paste0(FRM[1], paste0(FRM[3])))
          SPLIT_n <- strsplit(FRM2, ", n = (.*?),")[[1]] # change the length of iid random effects as data is different
          # recover length of each iid random effect groups
          nre_pr <- NULL
          if(is_Long & is.null(object[["REstrucS"]])){
            if(object$corLong) rmvCL = 0 else rmvCL=1
            for(nre_p in 1:length(object$famLongi)){
              if(nre_p<10){
                nre_10p = 0
              }else if(nre_p>=10){
                nre_10p = 1
              }
              nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
            }
            if(object$corLong) nre_pr <- sum(nre_pr)
          }
          if(object$corLong){
            FRM3 <- paste(paste(sapply(1:(1+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], ", n = ", nre_pr[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
          }else{
            FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], ", n = ", nre_pr[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
          }
          call.new <- object$call
          call.new[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
                                                           start=1,
                                                           stop=nchar(object$call[[length(object$call)]])-1),
                                                    ", dataOnly=TRUE, longOnly=TRUE)")
        }
        uData <- eval(parse(text=call.new)) # updated data with INLAjoint format
        if(!("E..coxph" %in% names(uData))){
          uData$E..coxph = rep(1, length(uData[[1]]))
        }
        # OTCrm <- sapply(object$longOutcome, function(x) grep(x, names(uData))) # exclude outfcomes
        # uData[-OTCrm] <- sapply(uData[-OTCrm], function(x) replace(x, is.na(x), 0), simplify=F)
        # rm(OTCrm)
        # if(strategy=="generic"){
        #   uData[-which(names(uData)==("Yjoint"))] <- sapply(uData[-which(names(uData)==("Yjoint"))], function(x) replace(x, is.na(x), 0), simplify=F)
        # }
        if(!is.list(uData)) uData <- as.list(as.data.frame(uData))
        nL_K <- length(uData[[1]])
        # now we prepare the precision matrix for all samples (large block diagonal matrix)
        IDshift <- 0
        if(exists("REnames") & strategy=="generic"){
          LengthUniqueID <- length(unique(na.omit(object$.args$data[[REnames[1]]])))
          uData[REnames] <- uData[paste0("W", object[["REstruc"]])]
        }
        if(!is.null(object[["REstrucS"]]) & strategy=="generic"){
          if(idLoop) uData[object[["REstrucS"]]] <- uData[paste0("W", substr(object[["REstrucS"]], 3, nchar(object[["REstrucS"]])))]
        }
        # A matrix for offset computation
        # ids to select the elements to keep in latent part of samples
        # baseline => substr(ct$tag, 1, 8)=="baseline" |
        A_off <- new("dgTMatrix", Dim=c(nL_K, sum(ct$length)))
        if(is_Long) A_off[, ct$start[SMPsel]] <- do.call(cbind, sapply(uData[ct$tag[SMPsel]], function(x) replace(x, is.na(x), 0), simplify=F))
        if(!is.null(object[["REstrucS"]]) | (is_Surv & strategy %in% c("dense", "split", "full", "full2", "full3", "joint"))){
          SMPselS <- which(ct$length==1 &
                             substr(ct$tag, nchar(ct$tag)-2, nchar(ct$tag)-1)=="_S" |
                             substr(ct$tag, nchar(ct$tag)-3, nchar(ct$tag)-2)=="_S") # if >10 markers
          if(length(SMPselS)>0){
            A_off[, ct$start[SMPselS]] <- do.call(cbind, sapply(uData[ct$tag[SMPselS]], function(x) replace(x, is.na(x), 0), simplify=F))
          }
        }
        # set baseline in A_off
        if(is_Surv & strategy!="joint"){
          Aproj <- NULL
          for(m in 1:M){
            if(object$basRisk[[m]] %in% c("rw1", "rw2")){
              colSEL <- ct$start[ct$tag==paste0("baseline", m, ".hazard")]:(ct$start[ct$tag==paste0("baseline", m, ".hazard")]+ct$length[ct$tag==paste0("baseline", m, ".hazard")]-1)
              linSEL <- which(!is.na(uData[[paste0("baseline", m, ".hazard.idx")]]))
              if(length(colSEL[na.omit(uData[[paste0("baseline", m, ".hazard.idx")]])])==1 & length(linSEL)==1){
                A_off[linSEL, colSEL[na.omit(uData[[paste0("baseline", m, ".hazard.idx")]])]] <- 1#diag(1, nrow=nrow(A_off[linSEL, colSEL[na.omit(uData[[paste0("baseline", m, ".hazard.idx")]])]]))
              }else{
                A_off[cbind(linSEL, colSEL[na.omit(uData[[paste0("baseline", m, ".hazard.idx")]])])] = 1#diag(1, nrow=nrow(A_off[linSEL, colSEL[na.omit(uData[[paste0("baseline", m, ".hazard.idx")]])]]))
              }
            }
          }
        }        # this is the offset used to evaluate the random effects
        offS <- (A_off %*% ParValMode)@x # mode (always use mode as first sample)
        if(strategy=="joint") assocNs <- NULL
        if(strategy %in% c("dense", "split", "full", "full2", "full3", "joint") & !is.null(assocNs)){
          if(NsampleREint==1){
            for(a_id in 1:length(assocNa)){
              # grab values of linear combination of fixed effects to share
              LP_sh <- offS[which(!is.na(uData$Yjoint[[grep(assocNa[a_id], names(uData$Yjoint))]]))]
              # scale it
              LP_shsc <- LP_sh * object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"]
              # add it to offset
              LPS_index <- which(!is.na(uData[[grep(paste0("^", assocNs[a_id], "$"), names(uData))]]))
              offS[LPS_index] <- offS[LPS_index] + LP_shsc
            }
          }else{
            offSet <- A_off %*% ParVal[, 1:(NsampleREint)] # samples
            for(a_id in 1:length(assocNa)){
              # grab values of linear combination of fixed effects to share
              LP_sh <- offSet[which(!is.na(uData$Yjoint[[grep(assocNa[a_id], names(uData$Yjoint))]])),]
              # scale it
              # LP_shsc <- LP_sh * c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
              #                      SMPH[1:(NsampleREint-1), grep(assocNs[a_id], colnames(SMPH))])
              if(!is.null(dim(LP_sh))){
                LP_shsc <- sapply(1:length(c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                             SMPH[1:(NsampleREint-1), grep(assocNs[a_id], colnames(SMPH))])),
                                  function(x) LP_sh[, x]*c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                                           SMPH[1:(NsampleREint-1), grep(assocNs[a_id], colnames(SMPH))])[x])
              }else{
                LP_shsc <- sapply(1:length(c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                             SMPH[1:(NsampleREint-1), grep(assocNs[a_id], colnames(SMPH))])),
                                  function(x) LP_sh[x]*c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                                           SMPH[1:(NsampleREint-1), grep(assocNs[a_id], colnames(SMPH))])[x])
              }
              # add it to offset
              LPS_index <- which(!is.na(uData[[grep(paste0("^", assocNs[a_id], "$"), names(uData))]]))
              offSet[LPS_index,] <- offSet[LPS_index, ] + LP_shsc
            }
          }
        }else{
          offSet <- A_off %*% ParVal[, 1:(NsampleREint)]
        }
        #Cmatrix <- as.matrix(BD_Cmat) # can we use sparse matrix here??
        #Cmatrix <- Matrix(BD_Cmat, sparse=T) # can we use sparse matrix here??
        # set up fixed residual errors for gaussian or lognormal families
        if(is_Long & strategy %in% c("generic", "dense", "split", "full", "full2", "full3", "joint")){
          ResErrScale <- rep(1, Nsample*nL_K)
          k_i <- 1
          if(is.null(object[["REstrucS"]])){
            for(k in 1:K){
              if(object$famLongi[k] %in% c("gaussian", "lognormal")){
                nL_k <- dim(ND)[1]
                posPrec <- which(substr(colnames(SMPH), 1, 18)=="Precision for the ")
                ResErrScale[rep(((k_i-1)*nL_k + 1):((k_i-1)*nL_k + nL_k), Nsample)+rep(nL_k*K*((1:Nsample)-1), each=nL_k)] <- rep(SMPH[, posPrec[k_i]], each=nL_k)
                k_i <- k_i+1
                ResErrFixed[[k]] <- list(hyper=list(prec=list(initial=0, fixed=TRUE)))
              }else{
                ResErrFixed[[k]] <- list()
              }
            }
          }else{
            ResErrFixed <- object$.args$control.family
          }
        }
        if(strategy=="generic" & is.null(object[["REstrucS"]])){
          # prepare A matrix (weights of the random effects)
          REweights <- INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))
          # when only slope random effect is included, it creats only zero rows in weight matrix
          # this is not accepted by inla call, therefore I add a tiny value at time zero for this
          # case to avoid the issue
          ZeroRE <- which(apply((INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))), 1, sum)==0)
          if(length(ZeroRE)>0 & is.null(object[["REstrucS"]]))  REweights[ZeroRE, which(REweights[ZeroRE+1,]!=0)] <- 1e-10
          A <- INLA::inla.as.sparse(kronecker(INLA::inla.as.sparse(diag(Nsample)), REweights))
          # A <- Diagonal(Nsample) %x% INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))
          # prepare outcome
          ncol_YJ <- ifelse(!is.null(object[["REstrucS"]]), length(uData$Yjoint), K)
          Yjoint <- matrix(unlist(uData$Yjoint), ncol=ncol_YJ)[rep(1:nL_K, Nsample),]
          IDnew <- 1:(nRE*Nsample)
          # fit the model to get random effects values for all samples
          if(NsampleRE!=F) SEL <- list("IDnew"=1:(nRE*Nsample)) else SEL<- NULL
          # NTH <- ifelse(parallel, "1:1", detectCores?)
          if(!is.null(object[["REstrucS"]])) FAM <- object$.args$family else FAM <- object$famLongi
          RE_estim <- INLA::inla(Yjoint ~ -1 + f(IDnew, model="generic0", Cmatrix=BD_Cmat,
                                                 hyper=list(theta=list(initial=0, fixed=TRUE))),
                                 family=FAM,
                                 data=list(Yjoint=Yjoint, IDnew=IDnew, A=A, offS=offS),
                                 offset=offS, scale=ResErrScale, #selection = SEL,
                                 control.predictor=list(A=A, link=1),
                                 control.compute=list(config=T),
                                 control.inla=list(int.strategy="eb"),
                                 control.family=ResErrFixed)
          # should I use empirical bayes here? is it worth?
          smpRE <- INLA::inla.posterior.sample(NsampleRE, RE_estim)
          RE_values <- matrix(unlist(sapply(smpRE, function(x) matrix(tail(x$latent, nRE*Nsample), nrow=nRE), simplify=F)), nrow=nRE)
          rm("smpRE")

        }else if(strategy %in% c("dense", "split", "full", "full2", "full3", "joint") & idLoopSet){ # full inla call
          RMNk <- object$.args$control.fixed$remove.names
          object$.args$control.fixed$remove.names <- c(object$.args$control.fixed$remove.names, rownames(object$summary.fixed))
          SEL <- NULL
          if(!is.null(nre_pr)){
            for(re_i in 1:length(nre_pr)){
              if(!is.null(object[["REstrucS"]])){
                if(re_i<=length(object[["REstrucS"]])){
                  SEL <- append(SEL, list((1:length(unique(ND[,id])))))
                }else{
                  for(re_j in 1:length(object[["REstruc"]])){
                    SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_j-1)))
                  }
                }
              }else{
                for(re_j in 1:length(grep(paste0("_L", re_i), object[["REstruc"]]))){
                  # for(re_j in 1:nre_pr[re_i]){
                  SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_j-1)))
                }
              }
            }
          }
          if(!is.null(object[["REstruc"]])) names_reL <-paste0("ID", object[["REstruc"]]) else names_reL <- NULL
          if(!is.null(object[["REstrucS"]])) names_reS <- object[["REstrucS"]] else names_reS <- NULL
          names(SEL) <- c(names_reS, names_reL)
          infoBHSEL <- 0
          if(is_Surv & strategy=="joint"){
            if(object$basRisk[[m]] %in% c("rw1", "rw2")){
              for(m in 1:M){
                SEL <- append(SEL, list(1:length(object$.args$data$baseline1.hazard.values)))
                infoBHSEL <- infoBHSEL + length(object$.args$data$baseline1.hazard.values)
              }
              names(SEL)[length(SEL)] <- paste0("baseline", m, ".hazard")
            }
          }
          if(NsampleREint>1) RE_values <- NULL
          if(is_Surv) BH_values <- NULL
          if(!SMPloop){ # unique inla() call for all samples and all individuals
            if(strategy=="full" | strategy=="full2" | strategy=="full3" | strategy=="joint"){
              if(NsampleREint==1){
                if(is_Surv & strategy!="joint"){
                  TETA <- data.frame(object$misc$theta.mode[-grep("baseline", object$misc$theta.tags)])
                }else{
                  TETA <- data.frame(object$misc$theta.mode)
                }
                # OFFSET <- data.frame(offS)
                uData <- append(uData, list("off"=as.matrix(offS)))
              }else{
                if(TRUE %in% sapply(c("rw1", "rw2"), function(x) x %in% unlist(object$basRisk)) & is_Surv & strategy!="joint"){
                  TETA <- sapply(1:NsampleREint, function(S) sapply(1:length(unname(SMPH[S,])), function(x) object$misc$to.theta[[x]](unname(SMPH[S,])[x]))[-grep("baseline", object$misc$theta.tags)])
                  if(is.null(dim(TETA))) TETA <- t(data.frame(TETA))
                }else{
                  TETA <- sapply(1:NsampleREint, function(S) sapply(1:length(unname(SMPH[S,])), function(x) object$misc$to.theta[[x]](unname(SMPH[S,])[x])))
                  if(is.null(dim(TETA))) TETA <- t(data.frame(TETA))
                }
                uData <- append(uData, list("off"=as.matrix(offSet)))
              }

              if(strategy=="full3"){
                TETA <- TETA[, seq(1, NsampleHY*NsampleFE, NsampleFE)]
                offS_NEW <- matrix(NA, nrow = length(uData[[1]]), ncol=NsampleHY)
                for(n_HY in 1:NsampleHY){
                  offS_HY <- uData$off[ ,1:NsampleFE + rep((n_HY-1)*NsampleFE, NsampleFE)]
                  for(m in 1:M){
                    # set survival part of offset for each Hyperpar sample (NsampleHY)
                    # and each FE sample (NsampleFE)
                    BH_m <- which(!is.na(uData[[paste0("expand", m, "..coxph")]]))
                    n_FEi <- length(BH_m)/NsampleFE
                    offS_NEW[BH_m, n_HY] <- c(sapply(1:NsampleFE, function(x) offS_HY[BH_m[1:n_FEi], x]))
                  }
                  LP_K <- sapply(object$longOutcome, function(x) grep(x, names(uData$Yjoint)))
                  for(k in 1:K){
                    # set longitudinal part of offset for each Hyperpar sample (NsampleHY)
                    # and each FE sample (NsampleFE)
                    if(length(LP_K[[k]])==0){
                      LP_k <- which(!is.na(uData$Yjoint))
                    }else{
                      LP_k <- which(!is.na(uData$Yjoint[[LP_K[k]]]))
                    }
                    n_FEi2 <- length(LP_k)/NsampleFE
                    offS_NEW[LP_k, n_HY] <- c(sapply(1:NsampleFE, function(x) offS_HY[LP_k[1:n_FEi2], x]))
                  }
                  AS_N <- sapply(assocNa, function(x) grep(x, names(uData$Yjoint)))
                  if(length(AS_N)>0){
                    for(n_asso in 1:length(AS_N)){
                      # set association part of offset for each Hyperpar sample (NsampleHY)
                      # and each FE sample (NsampleFE)
                      if(length(AS_N)==0){
                        AS_n <- which(!is.na(uData$Yjoint))
                      }else{
                        AS_n <- which(!is.na(uData$Yjoint[[AS_N[n_asso]]]))
                      }
                      n_FEi3 <- length(AS_n)/NsampleFE
                      offS_NEW[AS_n, n_HY] <- c(sapply(1:NsampleFE, function(x) offS_HY[AS_n[1:n_FEi3], x]))
                    }
                  }
                }
                uData$off <- offS_NEW
              }
              # INLA:::inla.tempdir()
              inla.setOption(malloc.lib='compiler')
              inla.setOption(INLAjoint.features=TRUE)
              object$.args$control.inla$compute.initial.values=FALSE
              wd <- INLA:::inla.tempdir()#"model.files"
              # unlink(wd, recursive = TRUE)
              if(length(which(object$.args$control.predictor$link!=1))>0) warning("Link function is not default, this has to be added here and has not yet been done. Please contact INLAjoint@gmail.com")
              if(!silentMode) message("Estimate conditional posterior of random effects...")
              if(Nsample==1) TETA <- TETA[,1]
              r <- inla(formula = formula(FRM3),
                        data = uData,
                        offset = off,
                        verbose = !TRUE,
                        working.directory = wd,
                        family = object$.args$family,
                        control.predictor=list(link=1), # to avoid warning when outcome includes NA
                        control.family = object$.args$control.family,
                        control.fixed = object$.args$control.fixed,
                        control.mode = list(theta = TETA,
                                            # x = object$mode$x,
                                            fixed = TRUE,
                                            restart = FALSE),
                        control.compute = list(return.marginals = FALSE, config=T),
                        control.inla = object$.args$control.inla,
                        quantiles = c(),
                        inla.call = "",
                        keep = TRUE,
                        safe = FALSE)
              if(strategy=="full3"){
                r <- INLA:::inla.run.many(NsampleHY, wd, num.threads = object$.args$num.threads, cleanup = !TRUE, verbose = !TRUE)#
              }else{
                r <- INLA:::inla.run.many(Nsample, wd, num.threads = object$.args$num.threads, cleanup = !TRUE, verbose = !TRUE)#
              }
              inla.setOption(INLAjoint.features=FALSE)
              inla.setOption(malloc.lib='mi')
              if(nRE==1 & length(unique(newData[, object$id]))==1 & strategy!="joint" & strategy!="full3"){
                RE_values <- unlist(sapply(1:Nsample, function(S) sapply(INLA::inla.posterior.sample(NsampleRE, r[[S]], selection=SEL), function(x) x$latent), simplify=F))
              }else if(strategy=="full3"){
                RE_values <- do.call(cbind, sapply(1:NsampleHY, function(S) sapply(INLA::inla.posterior.sample(NsampleRE, r[[S]], selection=SEL), function(x) x$latent), simplify=F))
                NRE_i <- length(SEL) # number of random effects
                NRE_ii <- (dim(RE_values)[1]/NRE_i)/NsampleFE # number of individuals
                id_REV <- sapply(1:NsampleFE, function(x) rep(1:NRE_ii, NRE_i) + rep((0:(NRE_i-1))*(NRE_ii*NsampleFE), each=NRE_ii) + rep(NRE_ii*(x-1), (NRE_ii*NRE_i)))
                RE_values <- do.call(cbind, apply(RE_values, 2, function(x) apply(id_REV, 2, function(xx) x[xx]), simplify=F))
              }else if(strategy!="joint"){
                RE_values <- do.call(cbind, sapply(1:Nsample, function(S) sapply(INLA::inla.posterior.sample(NsampleRE, r[[S]], selection=SEL), function(x) x$latent), simplify=F))
              }else if(strategy=="joint"){
                NAMs <- rownames(INLA::inla.posterior.sample(1, r[[1]], selection=SEL)[[1]]$latent)
                s_REV <- sapply(1:Nsample, function(S) sapply(INLA::inla.posterior.sample(NsampleRE, r[[S]], selection=SEL), function(x) x$latent), simplify=F)
                if(nRE==1 & length(unique(newData[, object$id]))==1){
                  RE_values <- c(sapply(s_REV, function(x) x[-grep("baseline", NAMs),]))#c((dim(s_REV[[1]])[1]-infoBHSEL+1):dim(s_REV[[1]])[1])
                }else if(nRE==1 & length(unique(newData[, object$id]))>1){
                }else if(nRE>1 & length(unique(newData[, object$id]))==1){
                }else if(nRE>1 & length(unique(newData[, object$id]))>1){
                  RE_values <- c(sapply(s_REV, function(x) x[-grep("baseline", NAMs),]))#c((dim(s_REV[[1]])[1]-infoBHSEL+1):dim(s_REV[[1]])[1])
                }else if(nRE==0){ # only baseline? or baseline + frailty? Treated as only baseline for now.

                }
                if(TRUE %in% sapply(object$basRisk, function(x) x%in%c("rw1", "rw2"))){
                  s_bas <- do.call(cbind, sapply(s_REV, function(x) x[c((dim(s_REV[[1]])[1]-infoBHSEL+1):dim(s_REV[[1]])[1]),], simplify=F))
                }
              }
            }else{
              for(REint in 1:NsampleREint){
                if(NsampleREint==1){
                  if(is_Surv){
                    TETA <- object$misc$theta.mode[-grep("baseline", object$misc$theta.tags)]
                  }else{
                    TETA <- object$misc$theta.mode
                  }
                  uData$offS <- offS
                }else{
                  if(TRUE %in% sapply(c("rw1", "rw2"), function(x) x %in% unlist(object$basRisk)) & is_Surv){
                    TETA <- sapply(1:length(unname(SMPH[REint,])), function(x) object$misc$to.theta[[x]](unname(SMPH[REint,])[x]))[-grep("baseline", object$misc$theta.tags)]
                  }else{
                    TETA <- sapply(1:length(unname(SMPH[REint,])), function(x) object$misc$to.theta[[x]](unname(SMPH[REint,])[x]))
                  }
                  uData$offS <- offSet[, REint]
                }
                RE_estim <- INLA::inla(formula = formula(FRM3),
                                       family = object$.args$family,
                                       data = uData,
                                       E = E..coxph,
                                       offset = offS,
                                       control.mode=list(theta=TETA, fixed=TRUE),
                                       control.family = object$.args$control.family,
                                       control.inla = object$.args$control.inla,
                                       control.fixed = object$.args$control.fixed,
                                       control.compute = list(config=T), verbose=F)
                if(strategy=="dense"){
                  # need to select baseline and use it for each sample instead of the ParVal one
                  smpRE <- sapply(INLA::inla.posterior.sample(NsampleRE, RE_estim, selection=SEL), function(x) x$latent)
                  if(NsampleREint==1){
                    if(is.null(dim(smpRE))){ # vector
                      RE_values <- smpRE#[rep(1:length(SEL), Nsample)]
                    }else{ # matrix
                      RE_values <- smpRE#[rep(1:length(SEL), Nsample),]
                    }
                  }else{
                    if(is.null(dim(smpRE))){ # vector
                      RE_values <- c(RE_values, smpRE)
                    }else{ # matrix
                      RE_values <- cbind(RE_values, smpRE)
                    }
                  }
                  # add parameter N sample hyperpar and fixed used for RE (only a fraction of total?) => need to compute theta...
                  rm("smpRE")
                }else if(strategy=="split"){
                  sapply(1:length(SEL), function(x) RE_estim$summary.random[[names(SEL)[x]]][SEL[[x]],"sd"])

                }
              }
            }

          }else{ # loop inla call over each sample (useful if this function crashes because the inla call is out of scale
            for(REint in 1:NsampleREint){
              if(NsampleREint==1){
                if(is_Surv){
                  TETA <- object$misc$theta.mode[-grep("baseline", object$misc$theta.tags)]
                }else{
                  TETA <- object$misc$theta.mode
                }
                uData$offS <- offS
              }else{
                if(TRUE %in% sapply(c("rw1", "rw2"), function(x) x %in% unlist(object$basRisk)) & is_Surv){
                  TETA <- sapply(1:length(unname(SMPH[REint,])), function(x) object$misc$to.theta[[x]](unname(SMPH[REint,])[x]))[-grep("baseline", object$misc$theta.tags)]
                }else{
                  TETA <- sapply(1:length(unname(SMPH[REint,])), function(x) object$misc$to.theta[[x]](unname(SMPH[REint,])[x]))
                }
                uData$offS <- offSet[, REint]
              }
              RE_estim <- INLA::inla(formula = formula(FRM3),
                                     family = object$.args$family,
                                     data = uData,
                                     E = E..coxph,
                                     offset = offS,
                                     control.mode=list(theta=TETA, fixed=TRUE),
                                     control.family = object$.args$control.family,
                                     control.inla = object$.args$control.inla,
                                     control.fixed = object$.args$control.fixed,
                                     control.compute = list(config=T), verbose=F)
              if(strategy=="dense"){
                # need to select baseline and use it for each sample instead of the ParVal one
                smpRE <- sapply(INLA::inla.posterior.sample(NsampleRE, RE_estim, selection=SEL), function(x) x$latent)
                if(NsampleREint==1){
                  if(is.null(dim(smpRE))){ # vector
                    RE_values <- smpRE#[rep(1:length(SEL), Nsample)]
                  }else{ # matrix
                    RE_values <- smpRE#[rep(1:length(SEL), Nsample),]
                  }
                }else{
                  if(is.null(dim(smpRE))){ # vector
                    RE_values <- c(RE_values, smpRE)
                  }else{ # matrix
                    RE_values <- cbind(RE_values, smpRE)
                  }
                }
                # add parameter N sample hyperpar and fixed used for RE (only a fraction of total?) => need to compute theta...
                rm("smpRE")
              }else if(strategy=="split"){
                sapply(1:length(SEL), function(x) RE_estim$summary.random[[names(SEL)[x]]][SEL[[x]],"sd"])
              }
            }
          }
          if(idLoopSet){ # save all rando effects before selecting for each individuals
            RE_valuesG <- RE_values
            idLoopSet <- FALSE
          }
          if(NsampleREint > Nsample){
            stop("Argument NsampleREint should be less or equal to Nsample")
          }else if(NsampleREint < Nsample){
            Nreps <- trunc(Nsample/NsampleREint)
            Nadds <- (Nsample %% NsampleREint)/NsampleREint
            if(is.null(dim(RE_values))){
              if(Nadds==0){
                AD_re <- NULL
              }else if(Nadds>0){
                AD_re <- trunc(1:(length(RE_values)*Nadds))
              }
              RE_values <- RE_values[c(rep(1:length(RE_values), Nreps), AD_re)]
              RE_valuesG <- RE_valuesG[c(rep(1:length(RE_valuesG), Nreps), AD_re)]
            }else{
              if(Nadds==0){
                AD_re <- NULL
              }else if(Nadds>0){
                AD_re <- trunc(1:(ncol(RE_values)*Nadds))
              }
              RE_values <- RE_values[, c(rep(1:ncol(RE_values), Nreps), AD_re)]
              RE_valuesG <- RE_valuesG[, c(rep(1:ncol(RE_valuesG), Nreps), AD_re)]
            }
          }
        }else if(strategy %in% c("dense", "split") & !is.null(object[["REstrucS"]])){
          # setup full model and fix hyperparameters
          # remove all fixed effects and use offset
          #estimate posterior random effects

          # RE_estim <- inla(Yjoint ~ )
        }
        if(is_Long){
          if(!idLoop){
            ND <- newData[newData[, object$id] == idPred,] # back to individuals now that random effects are done
            if(!is.null(dim(RE_valuesG))){ # only if there are more than 1 individual
              RE_values <- RE_valuesG[rep(1, length(object[["REstruc"]]))+length(unique(newData[,id]))*(seq(1, length(object[["REstruc"]]))-1)+(which(unique(newData[,id]) == idPred)-1),]
              RE_values <- RE_values[order(order(sapply(names(SEL), function(x) grep(x, ct$tag)))), ]
              # RE_values <- RE_valuesG[1:length(object[["REstruc"]])+c(1:length(object[["REstruc"]]))*which(unique(newData[,id]) == idPred),]
            }
          }
          if(FEonly) RE_values <- matrix(0, nrow = nrow(RE_values), ncol=ncol(RE_values))
          if(!is.null(object[["REstrucS"]])){
            RE_valuesL <- RE_values[-FRAIL_ind,]
          }else{
            RE_valuesL <- RE_values
          }
          # if(NsampleRE==F){ # when do we expect this to be FALSE?
          #   RE_values <- matrix(tail(RE_estim$summary.linear.predictor$mean, nRE*Nsample), nrow=nRE)
          # }
          # compute linear predictors for each sample at NtimePoints

          NEWdata[paste0("ID", object[["REstruc"]])] <- NEWdata[paste0("W", object[["REstruc"]])]
          # A matrix for offset computation
          A_LP <- new("dgTMatrix", Dim=c(length(NEWdata[[1]])-length(survPart), sum(ct$length)))
          if(K==1){
            Lout <- 1
          }else{
            Lout <- unique(c(sapply(object$longOutcome, function(x) grep(paste0(x, "_"), names(NEWdata$Yjoint)))))
          }
          indL <- unname(rep(1:NTP, length(Lout))+(rep(Lout, each=NTP)-1)*NTP)
          NEWdata <- suppressWarnings(sapply(NEWdata, function(x) replace(x, is.na(x), 0)))
          if(class(NEWdata)[1]=="matrix") NEWdata <- as.list(as.data.frame(NEWdata))
          if(!is.null(survPart)){
            A_LP[, ct$start[which(ct$tag %in% names(uData))]] <- sapply(ct$tag[which(ct$tag %in% names(uData))], function(x) NEWdata[[x]])[-survPart,]
          }else{
            A_LP[, ct$start[which(ct$tag %in% names(uData))]] <- sapply(ct$tag[which(ct$tag %in% names(uData))], function(x) NEWdata[[x]])
          }
          # paramVal with each configuration of the random effects
          if(NsampleRE!=F & loopRE){
            LP_long <- NULL
            for(NSRE in 1:NsampleRE){
              ParVal[posRE, ] <- RE_valuesL[, ((NSRE-1)*Nsample+1):((NSRE-1)*Nsample+Nsample)]
              LP_long <- rbind(LP_long, t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal)))
            }
          }else if(NsampleRE!=F){
            ParVal[posRE, ] <- 0
            ParVal2 <- ParVal[,rep(1:Nsample, each=NsampleRE)]
            if(!exists("RWBH")) RWBH <- 0
            if(is_Surv & strategy=="joint" & length(RWBH)>0){ # add baseline hazard random walk samples
              ParVal2[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                  (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                     ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- s_bas
            }
            POSc <- rep(1:Nsample, NsampleRE)+rep(0:(NsampleRE-1)*Nsample, each=Nsample)
            RE_mat <- new("dgTMatrix",
                          i = as.integer(rep(posRE, length(POSc))-1),
                          j = as.integer(rep(POSc, each=length(posRE))-1), x=c(RE_valuesL), Dim=dim(ParVal2))
            ParVal2 <- ParVal2 + RE_mat
            # ParVal2[posRE, POSc] <- RE_valuesL
            LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal2))
          }else{
            ParVal[posRE, ] <- RE_valuesL
            LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal))
          }
          # LP_long_t <- NULL
          if(F){
            errCT <- 1 # counter for error terms
            for(k in 1:K){
              if(!is.null(names(object$.args$data$Yjoint))){
                k_id <- grep(object$longOutcome[k], names(object$.args$data$Yjoint))[1]
              }else{
                k_id <- 1
              }
              LP_long_k <- LP_long[, (1:NTP)+(k_id-1)*NTP]
              long_i_den <- NULL
              if(!(NA %in% ND[, object$longOutcome[k]])){
                if(object$famLongi[k]=="gaussian"){

                  if(length(which(rep(TPO, K) %in% ND[, object$timeVar]))>1){
                    long_i_mu <- LP_long_k[, which(TPO %in% ND[, object$timeVar])]
                    long_i_true <- ND[, object$longOutcome[k]]
                    ResErr_i <- rep(sqrt(1/SMPH[, posPrec[errCT]]), each=NsampleRE)
                    long_i_den = c(long_i_den, sapply(1:dim(long_i_mu)[1],
                                                      function(c) prod(mapply(function(x,y) dnorm(x, mean=y, sd=ResErr_i[c]),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c,])))) # sum the logs
                    errCT <- errCT+1
                  }else{
                    long_i_mu <- LP_long_k
                    long_i_true <- ND[, object$longOutcome[k]]
                    ResErr_i <- rep(sqrt(1/SMPH[, posPrec[errCT]]), each=NsampleRE)
                    long_i_den = c(long_i_den, sapply(1:length(long_i_mu),
                                                      function(c) prod(mapply(function(x,y) dnorm(x, mean=y, sd=ResErr_i[c]),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c])))) # sum the logs
                    errCT <- errCT+1
                  }
                }else if(object$famLongi[k]=="poisson"){
                  if(length(which(rep(TPO, K) %in% ND[, object$timeVar]))>1){
                    long_i_mu <- LP_long_k[, which(TPO %in% ND[, object$timeVar])]
                    long_i_true <- ND[, object$longOutcome[k]]
                    long_i_den = c(long_i_den, sapply(1:dim(long_i_mu)[1],
                                                      function(c) prod(mapply(function(x,y) dpois(x, lambda=exp(y)),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c,])))) # sum the logs
                  }else{
                    long_i_mu <- LP_long_k
                    long_i_true <- ND[, object$longOutcome[k]]
                    long_i_den = c(long_i_den, sapply(1:length(long_i_mu),
                                                      function(c) prod(mapply(function(x,y) dpois(x, lambda=exp(y)),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c])))) # sum the logs
                  }
                }else if(object$famLongi[k]=="binomial"){
                  if(length(which(rep(TPO, K) %in% ND[, object$timeVar]))>1){
                    long_i_mu <- LP_long_k[, which(TPO %in% ND[, object$timeVar])]
                    long_i_true <- ND[, object$longOutcome[k]]
                    long_i_den = c(long_i_den, sapply(1:dim(long_i_mu)[1],
                                                      function(c) prod(mapply(function(x,y) dbinom(x, size=1, prob=exp(y)/(1+exp(y))),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c,])))) # sum the logs
                  }else{
                    long_i_mu <- LP_long_k
                    long_i_true <- ND[, object$longOutcome[k]]
                    long_i_den = c(long_i_den, sapply(1:length(long_i_mu),
                                                      function(c) prod(mapply(function(x,y) dbinom(x, size=1, prob=exp(y)/(1+exp(y))),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c])))) # sum the logs
                  }
                }
              }else{
                long_i_den <- rep(1, dim(LP_long_k)[1])
              }
              long_i_den3 <- c(sapply(1:Nsample, function(x) long_i_den[(x-1)*NsampleRE + 1:NsampleRE]/sum(long_i_den[(x-1)*NsampleRE + 1:NsampleRE])))
              # LP_long_save <- LP_long
              LP_long[, (1:NTP)+(k_id-1)*NTP] <- LP_long_k*long_i_den3
            }
            LP_long <- t(sapply(1:Nsample, function(x) colSums(LP_long[(x-1)*NsampleRE + 1:NsampleRE,])))
            LP_long <- LP_long[rep(1:dim(LP_long)[1], each=NsampleRE),] # this may not be a good idea...
          }
          if(resErrLong){ # add residual error
            famerr <- which(object$famLongi %in%c("gaussian", "lognormal"))
            hyperr <- sort(c(grep("Gaussian", colnames(SMPH)), grep("gaussian", colnames(SMPH)), grep("lognormal", colnames(SMPH))))
            if(length(famerr) != length(hyperr)) stop("There is an issue with residual error computations, please contact INLAjoint@gmail.com")
            REsamFull <- NULL
            for(fer in 1:length(famerr)){
              sdErr <- sqrt(1/c(object$summary.hyperpar[hyperr[fer], "0.5quant"], SMPH[-Nsample, hyperr[fer]])) # keep mode for first and then use samples
              REsam <- sapply(sdErr, function(x) rnorm(NsampleRE*NTP, mean=0, sd=x)) # residual error realizations
              # add residual errors to linear predictors
              REsamF <- do.call("rbind", lapply(1:ncol(REsam), function(x) matrix(REsam[,x], ncol = NTP)))
              LP_long[, indL] <- LP_long[, indL][, (1:NTP)+NTP*(famerr[fer]-1)] + REsamF
            }
          }
          if(return.samples){
            RESpredL <- t(LP_long[, indL])
            addNamesL <- paste0("Sample_", 1:ncol(RESpredL))
          }else{
            if(inv.link){
              namesLink <- object$misc$linkfunctions$names[Lout]
              LP_long2 <- LP_long
              for(k in 1:K){
                if(object$famLongi[k]=="lognormal"){
                  LP_long2[, indL][, 1:NTP + rep(NTP*(k-1), NTP)] <- INLA::inla.link.log(LP_long2[, indL][, 1:NTP + rep(NTP*(k-1), NTP)], inverse = TRUE)
                }else{
                  LP_long2[, indL][, 1:NTP + rep(NTP*(k-1), NTP)] <- eval(parse(text=paste0("inla.link.", namesLink[k])))(LP_long2[, indL][, 1:NTP + rep(NTP*(k-1), NTP)], inverse = TRUE)
                }
              }
              RESpredL <- t(apply(LP_long2[, indL], 2, SumStats))
              rm(LP_long2)
            }else{
              RESpredL <- t(apply(LP_long[, indL], 2, SumStats))
            }
            addNamesL <- c("Mean", "Sd", "quant0.025", "quant0.5", "quant0.975")
          }
          if(ChangeID){
            idPred2 <- newID[which(newID[,2]==idPred),1]
          }else{
            idPred2 <- idPred
          }
          newPredL <- data.frame(rep(rep(idPred2, length(LdataPred[, object$id])), K), rep(LdataPred[, object$timeVar], K),
                                 rep(object$longOutcome, each=NTP), RESpredL)
          colnames(newPredL) <- c(object$id, object$timeVar, "Outcome", addNamesL)
          predL <- rbind(predL, newPredL)
          rm("RE_valuesL")
        }
      }else if (strategy==4){

      }
    }
    ###          ###
    ### SURVIVAL ###
    ###          ###
    if(is_Surv){
      if(is_Long & is.null(Csurv)){
        CsurvSET <- max(ND[, object$timeVar])
      }else if(!is_Long & is.null(Csurv)){
        CsurvSET <- 0
      }
      startP <- ifelse(is.null(Csurv), CsurvSET, Csurv)  # start point for survival
      TPO2 <- TPO[TPO>=startP]
      NTP2 <- length(TPO2)
      NTP_s <- NTP-NTP2+1
      # survPart2 <- survPart[rep(which(TPO %in% TPO2), M)+rep(0:(M-1), each=NTP2)*NTP] # extract part where there is an actual risk
      survPart2 <- survPart[unlist(sapply(1:M, function(x) which(NEWdata[[paste0("baseline", x, ".hazard.time")]][NEWdata[[paste0("baseline", x, ".hazard.idx")]]!=0] %in% TPO2)))] # extract part where there is an actual risk
      # baseline risk setup
      if(baselineHaz=="PWconstant"){
        if(dim(ND)[1]==1){ # use existent cutpoints
          if(!is.null(object$dataSurv)) assign(paste0(object$dataSurv), ND)
          if(object$nameTimeSurv %in% colnames(ND)){
            ND[, object$nameTimeSurv] <- horizon
          }else{
            ND <- cbind(ND, horizon)
            colnames(ND)[length(colnames(ND))] <- object$nameTimeSurv
          }
          assign(paste0(object$dataSurv), ND)
          #Surv <- ND
          call.new3 <- object$call
          call.new3[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
                                                            start=1,
                                                            stop=nchar(object$call[[length(object$call)]])-1), TXT1,
                                                     ", dataOnly=TRUE)")
          NEWdata <- eval(parse(text=call.new3)) # maybe need to store functions of time in the object?
          ANewS <- matrix(0, ncol=length(paramVal), nrow=length(NEWdata$baseline1.hazard.idx))
          # FIX THE FOLLOWING FOR MULTIPLE SURVIVAL SUBMODELS??
          # here the diag should starts after the longitudinal, need to adjust this because it starts at the beginning for now
          diag(ANewS[, ct$start[ct$tag=="baseline1.hazard"]:(ct$start[ct$tag=="baseline1.hazard"]+length(NEWdata$baseline1.hazard.idx))]) <- 1
          ANewS[, ct$start[-which(ct$tag=="baseline1.hazard")]] <- do.call(cbind, NEWdata[ct$tag[-which(ct$tag=="baseline1.hazard")]])
          ANewS <- INLA::inla.as.sparse(ANewS, na.rm=TRUE)
          predSurv <- data.frame(ND[, object$id], NEWdata$baseline1.hazard.time, exp(matrix(c(as.matrix(ANewS %*% paramVal)), nrow=length(NEWdata$baseline1.hazard.idx), ncol=1)))
          colnames(predSurv) <- c(object$id, object$nameTimeSurv, paste0("LinPred_L", 1))
          return(predSurv)
        }else{ # compute risk value at new cutpoints

        }
      }else if(baselineHaz=="interpolation"){
        # for now we compute the risk from time zero, maybe don't need anything below Csurv? (bc no need for risk before "survival time")
        Aproj <- NULL
        for(m in 1:M){
          if(object$basRisk[[m]] %in% c("rw1", "rw2")){
            mesh1d <- INLA::inla.mesh.1d(loc = object$summary.random[[paste0("baseline", m, ".hazard")]]$ID, degree = 1)
            if(m==1){
              Aproj <- INLA::inla.spde.make.A(mesh = mesh1d, loc = NEWdata[[paste0("baseline", m, ".hazard.time")]][NEWdata[[paste0("baseline", m, ".hazard.idx")]]!=0][which(NEWdata[[paste0("baseline", m, ".hazard.time")]][NEWdata[[paste0("baseline", m, ".hazard.idx")]]!=0] %in% TPO2)])
            }else{
              Aproj <- bdiag(Aproj, INLA::inla.spde.make.A(mesh = mesh1d, loc = NEWdata[[paste0("baseline", m, ".hazard.time")]][which(NEWdata[[paste0("baseline", m, ".hazard.time")]][NEWdata[[paste0("baseline", m, ".hazard.idx")]]!=0] %in% TPO2)]))
            }
          }
        }
        # use weights in A to quantify uncertainty
      }else if(baselineHaz=="smooth"){
        Aproj <- diag(length(TPO2)*M)
      }
      if(is_Long){ # association
        # if(length(which(sapply(patternAsso, length)==0))>0) patternAsso <- patternAsso[-which(sapply(patternAsso, length)==0)] # exclude SRE_ind
        # need to set up the longitudinal shared part to scale with association parameters
        if(length(assocNs2)>0) LP_longs <- LP_long[, -indL][, rep(NTP_s:NTP, length(assocNs2))+rep(patternAsso2-1, each=NTP2)*NTP]
        # I assume all associations are contiguous here (I think it's always true!)
        ct2$start[assocPos] <- ct2$start[assocPos] - c(0, cumsum(ct2$length[assocPos])[-length(assocPos)]) + cumsum(c(0, rep(NTP2, length(assocNs)-1)))
        ct2$start[-c(1:assocPos[length(assocPos)])] <- ct2$start[-c(1:assocPos[length(assocPos)])] - sum(ct2$length[assocPos]) + NTP2*length(assocNs)#dim(LP_longs)[2]
        ct2$length[assocPos] <- NTP2
        # set association from NEWdata  to map ParamVal2 to ct2
        NEWdata[assocNs] <- sapply(NEWdata[assocNs], function(x) ifelse(x==0,0,1), simplify=F)
      }else{
        ct2 <- ct
      }
      A_SP <- new("dgTMatrix", Dim=c(length(survPart2), as.integer(sum(ct2$length))))
      # set up covariates for survival part (baseline and association done after)
      A_SP[, ct2$start[which(ct2$tag %in% names(NEWdata))]] <- sapply(ct2$tag[which(ct2$tag %in% names(NEWdata))], function(x) NEWdata[[x]][survPart2])
      # baseline
      BLpos <- which(ct2$tag%in%paste0("baseline", 1:M, ".hazard"))
      A_SP[, c(sapply(BLpos, function(x) ct2$start[x]:(ct2$start[x]+ct2$length[x]-1)))] <- Aproj
      if(is_Long){          # SET ASSOCIATION INDICATOR HERE INSTEAD OF IS_LONG
        # association
        # set up diagonal matrix for each association (corresponding to each time point)
        matAssoc <- A_SP[, ct2$start[assocPos][1]:(ct2$start[assocPos][1] + sum(ct2$length[assocPos])-1)]
        if(!is.null(dim(matAssoc))){
          assocPoints <- as.matrix(matAssoc[NTP2*(1:M-1)+1, seq(1, ncol(matAssoc), by=NTP2)])
        }else{
          assocPoints <- matAssoc
        }
        if(M==1){
          Addassoc <- do.call("cbind", sapply(1:length(assocPoints), function(x) Diagonal(NTP2), simplify=F))
        }else if(M>1){
          Addassoc <- kronecker(assocPoints, Diagonal(NTP2))
        }
        A_SP[, ct2$start[assocPos][1]:(ct2$start[assocPos][1] + sum(ct2$length[assocPos]) -1)] <- Addassoc
      }
      # merge
      # scale the association parameters
      if(NsampleRE==F) nsamplere <- 1 else nsamplere <- NsampleRE
      if(NsampleRE!=F & loopRE){
        LP_surv <- NULL
        for(NSRE in 1:NsampleRE){
          # need to fix this since the non loopRE version has been modified!
          # SET ASSOCIATION INDICATOR HERE INSTEAD OF IS_LONG
          if(is_Long) SASCP <- t(LP_longs[(1:Nsample + (NSRE-1)*Nsample), ] * sapply(assocNs, function(x) SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==x)])[, rep(1:length(assocPos), each=NTP2)])
          if(is_Long) ParValS <- rbind(ParVal[1:(ct$start[assocPos][1]-1), ], SASCP, ParVal[-c(1:(ct$start[assocPos][1] + sum(ct$length[assocPos]) -1)), ]) else ParValS <- ParVal
          LP_surv <- rbind(LP_surv, exp(t(as.matrix(INLA::inla.as.dgTMatrix(A_SP, na.rm=TRUE) %*% ParValS))))
        }
      }else{
        # set up matrix to have shared part from longitudinal (LP_longs) scaled by association parameters (assocScaler)
        if(is_Long){
          SASCP <- NULL
          DECAL <- 0 # need to shift when SRE_ind between two time dependent variables
          SRE_indic <- 1
          for(ias in 1:length(assocNs)){
            if(ias %in% SRE_inda){ # SRE_ind
              assocScaler <- SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==assocNs[ias])][rep(1:Nsample, each=nsamplere)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
              if(!is.null(dim(RE_values))){
                SASCP_t <- RE_values[SRE_indi[SRE_indic],]*assocScaler # time fixed so only 1 line required
                SRE_indic <- SRE_indic+1
              }else{
                SASCP_t <- RE_values*assocScaler # time fixed so only 1 line required
              }
              PZ <- which(ct2$tag==assocNs[ias]) # position of current assoc
              m_ind <- as.integer(strsplit(assocNs[ias], "_S")[[1]][2])
              # S_ind <- 1:NTP2+rep((m_ind-1), NTP2)*NTP2
              PRM <- (ct2$start[PZ]+1):(ct2$start[PZ]+(ct2$length[PZ]-1)) # remove other time points
              A_SP <- A_SP[, -PRM]
              A_SP[which(!is.na(A_SP[, ct2$start[PZ]])), ct2$start[PZ]] <- 1
              ct2$start[-c(1:PZ)] <- ct2$start[-c(1:PZ)] - length(PRM)
              ct2$length[PZ] <- 1
              DECAL <- DECAL + NTP2
            }else{ # other associations
              if(length(which(gsub("Beta for ", "", colnames(SMPH))==assocNs[ias]))>0){
                # assocScaler <- sapply(assocNs, function(x) SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==x)])[rep(1:Nsample, each=nsamplere), rep(1:length(assocPos), each=NTP2)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
                assocScaler <- SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==assocNs[ias])][rep(1:Nsample, each=nsamplere)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
                if(!is.null(dim(LP_longs))){
                  SASCP_t <- t(LP_longs[, (1:NTP2)+(ias-1)*NTP2-DECAL] * assocScaler)
                }else{
                  SASCP_t <- LP_longs * assocScaler
                }
              }else if(length(which(gsub(" \\(scopy mean\\)", "", gsub("Beta0 for NL_", "", colnames(SMPH)))==assocNs[ias]))>0){
                nb <- length(grep(assocNs[ias], colnames(SMPH)))
                prop <- INLAjoint.scopy.define(nb)
                k_NL <- as.integer(strsplit(strsplit(assocNs[ias], "_L")[[1]][2], "_S")[[1]][1])
                if(length(grep("CV", assocNs[ias])>0)){
                  x_NLid <- grep(paste0("uv", k_NL), names(object$summary.random))
                }else if(length(grep("CS", assocNs[ias])>0)){
                  x_NLid <- grep(paste0("us", k_NL), names(object$summary.random))
                }else if(length(grep("SRE", assocNs[ias])>0)){
                  x_NLid <- grep(paste0("usre", k_NL), names(object$summary.random))
                } # CV_CS not done here
                xval <- object$summary.random[[x_NLid]]$mean
                xx.loc <- min(xval) + (max(xval)-min(xval)) * (0:(nb - 1))/(nb - 1)
                iterSMP <- 0 # keep track of RE samples
                SASCP_t <- NULL
                for(nsmp in 1:Nsample){
                  funNL <- splinefun(xx.loc, prop$W %*% SMPH[nsmp, grep(assocNs[ias], colnames(SMPH))], method = "natural")
                  SASCP_t <- cbind(SASCP_t, t(apply(LP_longs[(1:nsamplere)+(nsamplere*iterSMP), (1:NTP2)+(ias-1)*NTP2], 2, function(x) x*funNL(x))))
                  iterSMP <- iterSMP+1
                }
              }
            }
            SASCP <- rbind(SASCP, SASCP_t)
          }
          ParValS <- rbind(ParVal[1:(ct$start[assocPos][1]-1), rep(1:Nsample, each=nsamplere)], SASCP, ParVal[-c(1:(ct$start[assocPos][1] + sum(ct$length[assocPos]) -1)), rep(1:Nsample, each=nsamplere)])
          if(strategy=="joint" & (TRUE %in% sapply(object$basRisk, function(x) x %in% c("rw1", "rw2")))){
            ParValS[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                 (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                    ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- s_bas
          }else if(strategy=="full2" & (TRUE %in% sapply(object$basRisk, function(x) x %in% c("rw1", "rw2")))){
            ParValS[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                 (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                    ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- SMPF$samples[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                                                                                                                        (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                                                                                                                           ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),]
          }
        }else{
          if(strategy=="joint" & (TRUE %in% sapply(object$basRisk, function(x) x %in% c("rw1", "rw2")))){
            ParValS <- ParVal[, rep(1:Nsample, each=nsamplere)]
            ParValS[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                  (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                     ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- s_bas
          }else if(strategy=="full2" & (TRUE %in% sapply(object$basRisk, function(x) x %in% c("rw1", "rw2")))){
            ParValS <- ParVal[, rep(1:Nsample, each=nsamplere)]
            ParValS[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                 (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                    ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),] <- SMPF$samples[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                                                                                                                        (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                                                                                                                           ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))),]
          }else{
            ParValS <- ParVal
          }
          PS_nosetup <- TRUE
        }
        if(!is.null(object[["REstrucS"]])){ # frailty terms?
          if(exists("PS_nosetup")) ParValS <- ParVal[, rep(1:ncol(ParVal), NsampleRE)]
          if(!is.null(dim(RE_values))){
            ParValS[ct$start[sapply(object[["REstrucS"]], function(x) which(ct2$tag==x))],] <- RE_values[FRAIL_ind,]
          }else{
            ParValS[ct$start[sapply(object[["REstrucS"]], function(x) which(ct2$tag==x))],] <- RE_values
          }
          m_inti1 <- which(ct2$tag %in% object[["REstrucS"]])
          for(m_intin in m_inti1){ # remove other time points
            PRM <- (ct2$start[m_intin]+1):(ct2$start[m_intin]+(ct2$length[m_intin]-1))
            A_SP <- A_SP[, -PRM]
            ParValS <- ParValS[-PRM,]
            ct2$start[-c(1:m_intin)] <- ct2$start[-c(1:m_intin)] - length(PRM)
            ct2$length[m_intin] <- 1
          }
          # shared frailty
          if(length(unlist(sapply(object[["REstrucS"]], function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))>0){
            for(ias in 1:length(unlist(sapply(object[["REstrucS"]], function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))){
              m_inti <- unlist(sapply(object[["REstrucS"]], function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH)))))[ias]
              m_ind <- na.omit(sapply(sapply(paste0(object[["REstrucS"]], "_S"), function(x) strsplit(colnames(SMPH)[m_inti], x)[[1]][2]), function(x) as.integer(x)))
              # S_ind <- 1:NTP2+rep((m_ind-1), NTP2)*NTP2#which(TPO %in% TPO2)+rep((m_ind-1), NTP2)*NTP
              PRM <- (ct2$start[m_inti]+1):(ct2$start[m_inti]+(ct2$length[m_inti]-1)) # remove other time points
              A_SP <- A_SP[, -PRM]
              ParValS <- ParValS[-PRM,]
              A_SP[which(!is.na(A_SP[, ct2$start[m_inti]])), ct2$start[m_inti]] <- 1
              ct2$start[-c(1:m_inti)] <- ct2$start[-c(1:m_inti)] - length(PRM)
              ct2$length[m_inti] <- 1
              # compute scaled frailty term and insert in shared part
              assocScaler <- SMPH[, m_inti][rep(1:Nsample, nsamplere)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
              REval_ind <- which(c(object[["REstruc"]], object[["REstrucS"]]) == strsplit(ct2$tag[m_inti], paste0("_S", m_ind))[[1]])
              if(!is.null(dim(RE_values))){
                ParValS[ct2$start[m_inti], ] <- RE_values[REval_ind, ]*assocScaler
              }else{
                ParValS[ct2$start[m_inti], ] <- RE_values*assocScaler
              }
            }
          }
        }
        LP_surv <- exp(t(as.matrix(INLA::inla.as.dgTMatrix(A_SP, na.rm=TRUE) %*% ParValS)))
      }
      if(return.samples){
        Risk12 <- t(LP_surv)
        NCol <- ncol(Risk12)
        addNamesS <- paste0("Sample_", 1:NCol)
      }else{
        NCol <- 5
        Risk12 <- t(apply(LP_surv, 2, SumStats))
        addNamesS <- c("Haz_Mean", "Haz_Sd", "Haz_quant0.025", "Haz_quant0.5", "Haz_quant0.975")
      }
      Risk13 <- matrix(0, nrow=NTP-NTP2, ncol=NCol)
      Surv13 <- matrix(1, nrow=NTP-NTP2, ncol=5)
      Risk2 <- NULL
      for(m in 1:M){
        Risk2 <- rbind(Risk2, rbind(Risk13, Risk12[1:NTP2 + rep(max(NTP2), NTP2)*(m-1),]))
      }
      if(is.null(object$timeVar)) TimeVar <- "time" else TimeVar <- object$timeVar
      if(ChangeID){
        idPred2 <- newID[which(newID[,2]==idPred),1]
      }else{
        idPred2 <- idPred
      }
      newPredS <- data.frame(rep(idPred2, M*NTP), rep(TPO, M),
                             rep(paste0("S_", 1:M), each=NTP), Risk2)
      colnames(newPredS) <- c(object$id, TimeVar, "Outcome", addNamesS)
      if(survival){
        # compute survival curve
        # take middle of intervals
        #TPsurv <- c(0, LdataPred[, object$timeVar][-1]-diff(LdataPred[, object$timeVar])/2)
        SurvSamp <- NULL
        for(m in 1:M){
          # startP <- ifelse(is.null(Csurv), max(ND[, object$timeVar]), Csurv) # start point for survival
          # TPsurv2 <- TPO[TPO>=startP]
          # rmTP <- c(length(TPO)*(m-1) + c(which(TPO<startP)), rep(1:length(TPO), M-1)+(rep((1:M)[-m], each=length(TPO))-1)*length(TPO))
          rmTP <- c(rep(1:length(TPO2), M-1)+(rep((1:M)[-m], each=length(TPO2))-1)*length(TPO2))
          #rmTP <- which(!(NEWdata[[paste0("baseline", m, ".hazard.time")]][survPart] %in% TPO2))





          # select only the time points non related to m to remove and then compute the value of "diff(TPO2)" for this part (in case of 1->10, 1->15 need to avoid diff for the junction...)

          if(length(rmTP)>0){
            SurvSamp2 <- exp(-apply(LP_surv[,-rmTP], 1, function(x) cumsum(x*c(0, diff(TPO2)))))
            SurvSampAdd <- matrix(1, nrow=length(which(TPO<startP)), ncol=dim(SurvSamp2)[2])
          }else{
            SurvSamp2 <- exp(-apply(LP_surv, 1, function(x) cumsum(x*c(0, diff(TPO2)))))
            SurvSampAdd <- NULL
          }
          SurvSamp <- rbind(SurvSamp, SurvSampAdd, SurvSamp2)
        }
        if(dim(SurvSamp)[1] != dim(newPredS)[1]){
          addF <- matrix(1, ncol=5, nrow=dim(newPredS)[1]-dim(SurvSamp)[1])
          SurvSampF <- rbind(addF, t(apply(SurvSamp, 1, SumStats)))
        }else{
          SurvSampF <- t(apply(SurvSamp, 1, SumStats))
        }
        newPredS <- cbind(newPredS, SurvSampF)
        colnames(newPredS)[(length(colnames(newPredS))-4):length(colnames(newPredS))] <- c("Surv_Mean", "Surv_Sd", "Surv_quant0.025", "Surv_quant0.5", "Surv_quant0.975")
      }
      if(CIF){
        CIFSamp <- vector("list", M)
        CIF_Samp_ <- NULL
        for(m in 1:M){
          rmTP <- c(rep(1:length(TPO2), M-1)+(rep((1:M)[-m], each=length(TPO2))-1)*length(TPO2))
          CIFSamp2 <- apply(LP_surv[,-rmTP], 1, function(x) cumsum(x*c(0, diff(TPO2))))
          CIFSamp[[m]] <- rbind(CIFSamp[[m]], CIFSamp2)
        }
        # compute overall survival
        oSurv <- sapply(1:dim(CIFSamp[[1]])[2], function(x) exp(-rowSums(sapply(1:M, function(m) CIFSamp[[m]][,x]))))
        for(m in 1:M){
          rmTP <- c(rep(1:length(TPO2), M-1)+(rep((1:M)[-m], each=length(TPO2))-1)*length(TPO2))
          if(length(rmTP)>0){
            CIFSamp_2 <- sapply(1:dim(LP_surv)[1], function(x) cumsum(LP_surv[x, -rmTP]*oSurv[, x]*c(0, diff(TPO2))))
            CIFSampAdd <- matrix(0, nrow=length(which(TPO<startP)), ncol=dim(CIFSamp_2)[2])
          }else{
            CIFSamp_2 <- sapply(1:dim(LP_surv)[1], function(x) cumsum(LP_surv[x, ]*oSurv[, x]*c(0, diff(TPO2))))
            CIFSampAdd <- NULL
          }
          CIF_Samp_ <- rbind(CIF_Samp_, CIFSampAdd, CIFSamp_2)
          if(dim(CIF_Samp_)[1] != dim(newPredS)[1]){
            addF <- matrix(0, ncol=5, nrow=dim(newPredS)[1]-dim(CIF_Samp_)[1])
            CIFSampF <- rbind(addF, t(apply(CIF_Samp_, 1, SumStats)))
          }else{
            CIFSampF <- t(apply(CIF_Samp_, 1, SumStats))
          }
        }
        newPredS <- cbind(newPredS, CIFSampF)
        colnames(newPredS)[(length(colnames(newPredS))-4):length(colnames(newPredS))] <- c("CIF_Mean", "CIF_Sd", "CIF_quant0.025", "CIF_quant0.5", "CIF_quant0.975")
      }
      predS <- rbind(predS, newPredS)
    }
    #list("PredL"=newPredL, "PredS"=newPredS)
  }
  # names(PRED) <- paste0("ID", unique(newData[, object$id]))
  # return(PRED)
  if(!silentMode) message(paste0("...done!"))
  return(list("PredL"=predL, "PredS"=predS))
}








