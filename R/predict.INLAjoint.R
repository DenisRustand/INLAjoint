#' @export
# @import doFuture
# @import foreach
# @import doRNG
# @import randtoolbox
#' @importFrom Matrix bdiag Diagonal
#' @importFrom methods new

predict.INLAjoint <- function(object, newData, timePoints=NULL, NtimePoints=50, strategy=1,
                              NsampleRE=50, loopRE=FALSE, Nsample=300, oldData=NULL, id=NULL, Csurv=NULL,
                              horizon=NULL, baselineHaz="interpolation", Nthreads=NULL,
                              return.samples=FALSE, parallel=FALSE, survival=FALSE, CIF=FALSE,
                              inv.link=FALSE, ...){ # oldData is temporary
  # id is the id column name in dataset for survival data only (otherwise it's given by longitudinal)
  # Csurv is to get predictions conditional on survival up to given time
  # strategy: 1=default ; 2=full sampling ; 3=update ; 4=analytical
  if (!"INLAjoint" %in% class(object)){
    stop("Please provide an object of class 'INLAjoint' (obtained with joint() function).\n")
  }
  # baselineHaz = "smooth" | "interpolation"
  out <- NULL
  SumStats <- function(x) return(c(mean(x), sd(x), quantile(x, c(0.025, 0.5, 0.975))))
  if(!is.null(object$id)) id <- object$id else if(is.null(id)) stop("Please specify individual id column name with argument 'id'")
  is_Long <- is_Surv <- FALSE
  if(!is.null(object$REstruc)){
    idVect <- na.omit(unique(object$.args$data[[paste0("ID", object$REstruc[[1]])]]))
    is_Long <- TRUE
  }
  if(!is.null(object$SurvInfo)){
    if(!exists("idVect")){
      idVect <- unique(object$.args$data$expand1..coxph)
    } else{
      if(!any(idVect %in% unique(object$.args$data$expand1..coxph))) stop("id mismatch between longi and surv.")
    }
    is_Surv <- TRUE
    M <- length(object$survOutcome) # number of survival outcomes
    # check if we have baseline info for horizon
    for(m in 1:M){
      # add " method is interpolatiioon and forecast is required => switching to smooth
      # explain that interpolation means constant after last time point and suggest to use smooth to use a smooth prediction of baseline after
      if(horizon>max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]) & baselineHaz=="interpolation"){
        warning(paste0("The fitted model has baseline risk information up until value ",
                       max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]), " for survival outcome ", m, ". Since you ask for prediction at horizon ", horizon, " I will assume constant baseline hazard beyond the maximum available value. Alternatively, you can use baselineHaz='smooth' to use splines to predict the baseline hazard (for each sample). Alternatively, adding 'horizon' in the control options of the inla() call allows to extend the baseline beyond the last observed event time (linear extension, less flexible than the smooth method)."))
      }
    }
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
  if(is.null(timePoints)){
    if(is.null(horizon)){
      timePoints <- seq(0, max(newData[, object$timeVar]), len=NtimePoints)
    }else{#} if(Csurv==0){
      timePoints <- seq(0, horizon, len=NtimePoints)
      # }else{
      #need to have a time point at Csurv there
    }
  }
  firstID <- unique(newData[, object$id])[1]
  if(strategy %in% c(1,2)){
    message("Start sampling")
    SMPH <- INLA::inla.hyperpar.sample(Nsample, object)
    SMP <- INLA::inla.rjmarginal(Nsample, object)
    message("Sampling done.")
  }
  ParVal <- new("dgTMatrix", Dim=c(sum(ct$length), as.integer(Nsample)))
  if(is_Surv){
    ParVal[c(c(sapply(1:M, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                        (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                           ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
             ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                              !ct$tag %in% paste0("baseline", 1:M, ".hazard"))]),] <- SMP$samples#[which(!substr(rownames(SMP$samples),1, 9) %in% paste0("baseline", 1:M)),]
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
    ParVal[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))]),] <- SMP$samples
  }
  if(is_Long){
    K <- length(object$famLongi) # number of longitudinal outcomes
    ct2 <- ct
    call.new <- object$call
    call.new[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
                                                  start=1,
                                                  stop=nchar(object$call[[length(object$call)]])-1),
                                           ", dataOnly=TRUE, longOnly=TRUE)")
    lenPV <- length(paramVal)
    SMPsel <- which(ct$length==1 &
                      substr(ct$tag, nchar(ct$tag)-2, nchar(ct$tag)-1)=="_L" |
                      substr(ct$tag, nchar(ct$tag)-3, nchar(ct$tag)-2)=="_L") # if >10 markers
    NamesH <- colnames(SMPH)
    nRE <- length(object$REstruc)
    if(nRE==1){
      BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample))) # adapt size
      diag(BD_Cmat) <- sqrt(1/SMPH[, which(substr(colnames(SMPH), 1, 16)=="Precision for ID")])
    }else if(nRE>1){
      # identify the position of the cholesky elements in hyperparameters
      if(object$corLong){
        PosH <- which(substr(NamesH, 1, 5)=="Theta" &
                        substr(NamesH, nchar(NamesH)-nchar(object$REstruc[1])-1,
                               nchar(NamesH))==paste0("ID", object$REstruc[1]))
        if(strategy%in%c(1,2)){
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
          # browser()
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
                          substr(NamesH, nchar(NamesH)-nchar(object$REstruc[nRE_pk])-1,
                                 nchar(NamesH))==paste0("ID", object$REstruc[nRE_pk]))
          nRE_k <- length(which(substr(object$REstruc, nchar(object$REstruc)-2, nchar(object$REstruc))==paste0("_L", k) |
                                  substr(object$REstruc, nchar(object$REstruc)-3, nchar(object$REstruc))==paste0("_L", k)))
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
    # RE_values2 <- mvtnorm::rmvnorm(NsampleRE, sigma=solve(BD_Cmat))
    # RE_values <-t(sapply(1:nRE, function(x) RE_values2[, seq(x, Nsample*nRE, by=nRE)]))
    ResErrFixed <- vector("list", K)
    REnames <- paste0("ID", object$REstruc)
    posRE <- ct$start[sapply(REnames, function(x) which(ct$tag==x))]
    if(is_Surv){
      assocNs <- object$assoc
      assocPos <- sapply(assocNs, function(x) which(ct$tag==x))
      # identify the longitudinal needed for association
      # first identify shared part from longitudinal (no duplicates, so if CV from longitudinal 1 is shared twice, we need to repeat it)
      outcomeAssoc <- names(object$.args$data$Yjoint)[which(substr(names(object$.args$data$Yjoint), 1, nchar(names(object$.args$data$Yjoint))-1) %in% ct2$tag)]
      outcomeAssoc2 <- sapply(outcomeAssoc, function(x) strsplit(x, split = "_S")[[1]][1])
      requiredAssoc <- sapply(assocNs, function(x) strsplit(x, split = "_S")[[1]][1])
      patternAsso <- unname(sapply(requiredAssoc, function(x) which(x==outcomeAssoc2)))
    }
  }
  # if(length(unique(newData[, object$id]))==1) Nthreads=1
  # if(is.null(Nthreads)) Nthreads = future::availableCores("mc.cores")-1
  # registerDoFuture()
  # registerDoRNG()
  # message(paste0("Number of threads: ", Nthreads))
  # plan("multisession", workers = Nthreads)
  # PRED <- foreach(idPred = unique(newData[, object$id]), .export=c(ls(), ls(envir=.GlobalEnv), lsf.str(envir=.GlobalEnv), "joint")) %dopar% {
  for(idPred in unique(newData[, object$id])){
    message(paste0("Computing longitudinal predictions for individual ", idPred, " on PID: ", Sys.getpid()))
    ND <- newData[newData[, object$id] == idPred,] # select only one id at a time
    if(is_Long & strategy==2){ # add time points at observed longi to compute density
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
      horizon2 <- max(TPO)+0.0001
      SdataPred <- ND[!duplicated(ND[, object$id]),]
      #if(!is.null(object$dataSurv))
      for(m in 1:M){
        if(length(paste0(object$SurvInfo[[m]]$nameTimeSurv))>1) NTS <- tail(paste0(object$SurvInfo[[m]]$nameTimeSurv),1) else NTS <- paste0(object$SurvInfo[[m]]$nameTimeSurv)
        if(length(paste0(object$SurvInfo[[m]]$survOutcome))>1) SVO <- tail(paste0(object$SurvInfo[[m]]$survOutcome),1) else SVO <- paste0(object$SurvInfo[[m]]$survOutcome)
        if(NTS %in% colnames(SdataPred)){
          SdataPred[, NTS] <- horizon2
        }else{
          SdataPred <- cbind(SdataPred, horizon2)
          colnames(SdataPred)[length(colnames(SdataPred))] <- NTS
        }
        if(!(SVO %in% colnames(SdataPred))){
          SdataPred <- cbind(SdataPred, 0)
          colnames(SdataPred)[length(colnames(SdataPred))] <- SVO
        }
      }
      assign(paste0(object$dataSurv), SdataPred)
      CTP <- c(TPO, max(TPO)+tail(diff(TPO),1)) # add fake point to extend and have the last value
      if(!is.null(object$cutpoints)) assign(paste0(object$cutpoints), CTP) else TXT1 <- ", cutpoints=CTP"
    }
    LdataPred <- ND[rep(1, length(TPO)), ]
    LdataPred[, object$timeVar] <- TPO
    if(!is.null(object$dataLong)) assign(paste0(object$dataLong), LdataPred)
    if(!is.null(object$dataSurv) & is_Surv) assign(paste0(object$dataSurv), SdataPred)
    call.new2[[length(object$call)]] <- paste0(substr(object$call[[length(object$call)]],
                                                   start=1,
                                                   stop=nchar(object$call[[length(object$call)]])-1), TXT1,
                                            ", dataOnly=TRUE)")
    NEWdata <- suppressWarnings(eval(parse(text=call.new2))) # maybe need to store functions of time in the object?
    survPart <- NULL
    if(is_Surv) survPart <- c(sapply(1:M, function(x) which(!is.na(eval(parse(text=paste0("NEWdata$baseline", x, ".hazard.idx")))))))
    if(!is.list(NEWdata)) NEWdata <- as.list(as.data.frame(NEWdata))
    if(is_Long){
      ###              ###
      ### LONGITUDINAL ###
      ###              ###
      ctL <- ct
      if(idPred %in% idVect & identical(ND[, object$timeVar], object$.args$data[[paste0(object$timeVar, "_L1")]][which(object$.args$data$IDIntercept_L1==idPred)])){
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
        # if(!is.null(object$REstruc)){
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
        # NEWdata[REnames] <- NEWdata[paste0("W", object$REstruc)]
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

      }else if(strategy %in% c(1,2)){
        # convert data into INLAjoint format with dataOnly option
        assign(paste0(object$dataLong), ND)
        uData <- eval(parse(text=call.new)) # updated data with INLAjoint format
        # OTCrm <- sapply(object$longOutcome, function(x) grep(x, names(uData))) # exclude outfcomes
        # uData[-OTCrm] <- sapply(uData[-OTCrm], function(x) replace(x, is.na(x), 0), simplify=F)
        # rm(OTCrm)
        uData[-which(names(uData)==("Yjoint"))] <- sapply(uData[-which(names(uData)==("Yjoint"))], function(x) replace(x, is.na(x), 0), simplify=F)
        if(!is.list(uData)) uData <- as.list(as.data.frame(uData))
        nL_K <- length(uData[[1]])
        # A matrix for offset computation
        # ids to select the elements to keep in latent part of samples
        # baseline => substr(ct$tag, 1, 8)=="baseline" |
        A_off <- new("dgTMatrix", Dim=c(nL_K, sum(ct$length)))
        A_off[, ct$start[SMPsel]] <- do.call(cbind, uData[ct$tag[SMPsel]])
        offS <- (A_off %*% ParVal)@x # this is the offset used to evaluate the random effects
        # now we prepare the precision matrix for all samples (large block diagonal matrix)
        IDshift <- 0
        LengthUniqueID <- length(unique(na.omit(object$.args$data[[REnames[1]]])))
        uData[REnames] <- uData[paste0("W", object$REstruc)]

        #Cmatrix <- as.matrix(BD_Cmat) # can we use sparse matrix here??
        #Cmatrix <- Matrix(BD_Cmat, sparse=T) # can we use sparse matrix here??
        # set up fixed residual errors for gaussian or lognormal families
        ResErrScale <- rep(1, Nsample*nL_K)
        k_i <- 1
        for(k in 1:K){
          if(object$famLongi[k] %in% c("gaussian", "lognormal")){
            nL_k <- dim(ND)[1]
            posPrec <- which(substr(colnames(SMPH), 1, 18)=="Precision for the ")
            ResErrScale[rep(((k_i-1)*nL_k + 1):((k_i-1)*nL_k + nL_k), Nsample)+rep(nL_k*K*((1:Nsample)-1), each=nL_k)] <- sqrt(1/rep(SMPH[, posPrec[k_i]], each=nL_k))
            k_i <- k_i+1
            ResErrFixed[[k]] <- list(hyper=list(prec=list(initial=0, fixed=TRUE)))
          }else{
            ResErrFixed[[k]] <- list()
          }
        }
        if(strategy==1){
          # prepare A matrix (weights of the random effects)
          REweights <- INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))
          # when only slope random effect is included, it creats only zero rows in weight matrix
          # this is not accepted by inla call, therefore I add a tiny value at time zero for this
          # case to avoid the issue
          ZeroRE <- which(apply((INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))), 1, sum)==0)
          if(length(ZeroRE)>0)  REweights[ZeroRE, which(REweights[ZeroRE+1,]!=0)] <- 1e-10
          A <- INLA::inla.as.sparse(kronecker(INLA::inla.as.sparse(diag(Nsample)), REweights))
          # A <- Diagonal(Nsample) %x% INLA::inla.as.sparse(matrix(unlist(uData[REnames]), ncol=nRE))
          # prepare outcome
          Yjoint <- matrix(unlist(uData$Yjoint), ncol=K)[rep(1:nL_K, Nsample),]
          IDnew <- 1:(nRE*Nsample)
          # fit the model to get random effects values for all samples
          if(NsampleRE!=F) SEL <- list("IDnew"=1:(nRE*Nsample)) else SEL<- NULL
          # NTH <- ifelse(parallel, "1:1", detectCores?)
          RE_estim <- INLA::inla(Yjoint ~ -1 + f(IDnew, model="generic0", Cmatrix=BD_Cmat,
                                         hyper=list(theta=list(initial=0, fixed=TRUE))),
                         family=object$famLongi,
                         data=list(Yjoint=Yjoint, IDnew=IDnew, A=A, offS=offS),
                         offset=offS, scale=ResErrScale, #selection = SEL,
                         control.predictor=list(A=A, link=1),
                         control.compute=list(config=T),
                         control.inla=list(int.strategy="eb"),
                         control.family=ResErrFixed)#, num.threads="1:1")#, safe =F)
          # can I use empirical bayes here? is it worth?
          smpRE <- INLA::inla.posterior.sample(NsampleRE, RE_estim)
          RE_values <- matrix(unlist(sapply(smpRE, function(x) matrix(tail(x$latent, nRE*Nsample), nrow=nRE), simplify=F)), nrow=nRE)
          rm("smpRE")
        }
        # if(NsampleRE==F){ # when do we expect this to be FALSE?
        #   RE_values <- matrix(tail(RE_estim$summary.linear.predictor$mean, nRE*Nsample), nrow=nRE)
        # }
        # compute linear predictors for each sample at NtimePoints

        NEWdata[REnames] <- NEWdata[paste0("W", object$REstruc)]
        # A matrix for offset computation
        A_LP <- new("dgTMatrix", Dim=c(length(NEWdata[[1]])-length(survPart), sum(ct$length)))
        if(K==1){
          Lout <- 1
        }else{
          Lout <- sapply(object$longOutcome, function(x) grep(x, names(NEWdata$Yjoint)))
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
            ParVal[posRE, ] <- RE_values[, ((NSRE-1)*Nsample+1):((NSRE-1)*Nsample+Nsample)]
            LP_long <- rbind(LP_long, t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal)))
          }
        }else if(NsampleRE!=F){
          ParVal[posRE, ] <- 0
          ParVal2 <- ParVal[,rep(1:Nsample, NsampleRE)]
          POSc <- rep(1:Nsample, NsampleRE)+rep(0:(NsampleRE-1)*Nsample, each=Nsample)
          RE_mat <- new("dgTMatrix",
                    i = as.integer(rep(posRE, length(POSc))-1),
                    j = as.integer(rep(POSc, each=length(posRE))-1), x=c(RE_values), Dim=dim(ParVal2))
          ParVal2 <- ParVal2 + RE_mat
          # ParVal2[posRE, POSc] <- RE_values
          LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal2))
        }else{
          ParVal[posRE, ] <- RE_values
          LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal))
        }
        if(strategy==2){
          for(k in 1:K){
            long_i_den <- NULL
            if(object$famLongi[k]=="gaussian"){
              long_i_mu <- LP_long[, indL][, which(rep(TPO, K) %in% ND[, object$timeVar])][,(k-1)*nL_k + 1:nL_k]
              long_i_true <- ND[, object$longOutcome[k]]
              ResErr_i <- rep(sqrt(1/SMPH[, posPrec[k]]), each=NsampleRE)
              long_i_den = c(long_i_den, sapply(1:dim(long_i_mu)[1],
                                  function(c) prod(mapply(function(x,y) dnorm(x, mean=y, sd=ResErr_i[c]),
                                                                              x = long_i_true,
                                                                              y = long_i_mu[c,])))) # sum the logs
            }
          }
          LP_long_save <- LP_long
          A <- LP_long*long_i_den
          LP_long <- t(sapply(1:Nsample, function(x) colSums(A[(x-1)*NsampleRE + 1:NsampleRE,])))
          LP_long <- LP_long[rep(1:dim(LP_long)[1], each=NsampleRE),]
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
        newPredL <- data.frame(rep(LdataPred[, object$id], K), rep(LdataPred[, object$timeVar], K),
                               rep(object$longOutcome, each=NTP), RESpredL)
        colnames(newPredL) <- c(object$id, object$timeVar, "Outcome", addNamesL)
        # browser()

        predL <- rbind(predL, newPredL)
      }else if (strategy==4){

      }
    }
    ###          ###
    ### SURVIVAL ###
    ###          ###
    if(is_Surv){
      message(paste0("Computing survival predictions for individual ", idPred))
      if(is_Long & is.null(Csurv)){
        CsurvSET <- max(ND[, object$timeVar])
      }else if(!is_Long & is.null(Csurv)){
        CsurvSET <- 0
      }
      startP <- ifelse(is.null(Csurv), CsurvSET, Csurv)  # start point for survival
      TPO2 <- TPO[TPO>=startP]
      NTP2 <- length(TPO2)
      NTP_s <- NTP-NTP2+1
      survPart2 <- survPart[rep(which(TPO %in% TPO2), M)+rep(0:(M-1), each=NTP2)*NTP] # extract part where there is an actual risk
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
          mesh1d <- INLA::inla.mesh.1d(loc = object$summary.random[[paste0("baseline", m, ".hazard")]]$ID, degree = 1)
          if(m==1){
            Aproj <- INLA::inla.spde.make.A(mesh = mesh1d, loc = TPO2)
          }else{
            Aproj <- bdiag(Aproj, INLA::inla.spde.make.A(mesh = mesh1d, loc = TPO2))
          }
        }
        # use weights in A to quantify uncertainty
      }else if(baselineHaz=="smooth"){
        Aproj <- diag(length(TPO2)*M)
      }
      if(is_Long){ # association
        if(length(which(sapply(patternAsso, length)==0))>0) patternAsso <- patternAsso[-which(sapply(patternAsso, length)==0)] # exclude SRE_ind
        # need to set up the longitudinal shared part to scale with association parameters
        LP_longs <- LP_long[, -c(indL, survPart)][, rep(NTP_s:NTP, length(assocNs))+rep(patternAsso-1, each=NTP2)*NTP]
          # I assume all associations are contiguous here (I think it's always true!)
        ct2$start[assocPos] <- ct2$start[assocPos] - cumsum(c(0, ct2$length[assocPos][-1])) + cumsum(c(0, rep(NTP2, length(assocNs)-1)))
        ct2$start[-c(1:assocPos[length(assocPos)])] <- ct2$start[-c(1:assocPos[length(assocPos)])] - sum(ct2$length[assocPos]) + dim(LP_longs)[2]
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
        assocPoints <- as.matrix(matAssoc[NTP2*(1:M-1)+1, seq(1, ncol(matAssoc), by=NTP2)])
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
        if(is_Long) assocScaler <- sapply(assocNs, function(x) SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==x)])[rep(1:Nsample, nsamplere), rep(1:length(assocPos), each=NTP2)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
        if(is_Long) SASCP <- t(LP_longs * assocScaler)
        if(is_Long) ParValS <- rbind(ParVal[1:(ct$start[assocPos][1]-1), rep(1:Nsample, nsamplere)], SASCP, ParVal[-c(1:(ct$start[assocPos][1] + sum(ct$length[assocPos]) -1)), rep(1:Nsample, nsamplere)]) else ParValS <- ParVal
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
      newPredS <- data.frame(rep(idPred, M*NTP), rep(TPO, M),
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
  return(list("PredL"=predL, "PredS"=predS))
}








