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
#' @param startTime define a starting time for predictions.
#' @param horizon horizon of the prediction.
#' @param baselineHaz method used to evaluate the baseline hazard value, default is 'interpolation'
#' which is currently recommended. Experimental alternatives are being developed, including 'splines'
#' for an interpolation with splines but has not been properly validated with simulations yet.
#' @param return.samples boolean, when set to TRUE the samples are returned instead of summary
#' statistics over the samples. Default is FALSE.
#' @param FEonly boolean, when set to TRUE, only fixed effects are involved for predictions computations.
#' @param survival boolean, when set to TRUE the summary statistics over survival functions are
#' computed in addition to the summary statistics over the risk functions.
#' @param CIF boolean, when set to TRUE the summary statistics over cumulative incidence functions are
#' computed in addition to the summary statistics over the risk functions. Only applies to competing risks.
#' @param inv.link boolean, when set to TRUE the summary statistics are computed over the predictions of
#' longitudinal components after applying the inverse link function for each samples in addition to the
#' summary statistics over the linear predictors.
#' @param NidLoop Gives the number of individuals for which we compute predictions at once. For large number
#' of individuals, this will loop over groups of 'NidLoop' individuals and could make predictions computations faster.
#' @param resErrLong boolean, when set to TRUE the residual error for Gaussian or lognormal longitudinal
#' outcomes is added to the uncertainty of predictions (default is FALSE which predicts the true underlying
#' value of the longitudinal marker, i.e., error-free).
#' @param set.samples replace random effects with pre-sampled values.
#' @param silentMode a boolean that will stop printing messages during computations if turned to TRUE.
#' @param ... Extra arguments.
#' @export
#' @importFrom Matrix bdiag Diagonal
#' @importFrom methods new

predict.INLAjoint <- function(object, newData=NULL, newDataSurv=NULL, timePoints=NULL, NtimePoints=50,
                              NsampleHY=20, NsampleFE=20, NsampleRE=50, id=NULL, Csurv=NULL, startTime=NULL,
                              horizon=NULL, baselineHaz="interpolation", return.samples=FALSE, FEonly=FALSE,
                              survival=FALSE, CIF=FALSE, inv.link=FALSE, NidLoop="auto", resErrLong=FALSE,
                              set.samples=NULL, silentMode=FALSE, ...){
  # idGroup: loop over groups over random effects (useful if scaling issues)
  arguments <- list(...)
  # id is the id column name in dataset for survival data only (otherwise it's given by longitudinal)
  # Csurv is to get predictions conditional on survival up to given time
  idLoop=FALSE
  REmsg <- TRUE
  if(exists("object$run")) if(!object$run) stop("Please run the model (with function `joint.run()`)")
  if(is.null(newData)){ # if no new data is provided, return predicted fitted values
    PRED <- object$summary.fitted.values
    OUtc <- as.data.frame(object$.args$data$Yjoint)
    PRED$Outcome <- sapply(1:dim(PRED)[1], function(x) colnames(OUtc)[which(!is.na(OUtc[x,]))])
    return(PRED)
  }
  # Identify spatial hyperparameters (specific patterns to avoid false positives with regular IID effects)
  is_spatial_hyperpar <- function(name) {
    any(grepl("Phi for ID|Range for ID|Stdev for ID", name))
  }
  # Function to detect if a random effect is spatial
  is_spatial_re <- function(re_name, spatial_hyperpar_names) {
    if(is.null(re_name) || is.na(re_name) || re_name == "") return(FALSE)
    if(length(spatial_hyperpar_names) == 0) return(FALSE)
    if(any(is.na(spatial_hyperpar_names)) || any(spatial_hyperpar_names == "")) return(FALSE)
    tryCatch({
      any(sapply(spatial_hyperpar_names, function(s) {
        if(is.null(s) || is.na(s) || s == "") return(FALSE)
        grepl(s, re_name, fixed = TRUE)
      }))
    }, error = function(e) {
      return(FALSE)
    })
  }
  loopRE <- FALSE
  if(!"INLAjoint" %in% class(object)){
    stop("Please provide an object of class 'INLAjoint' (obtained with joint() function).\n")
  }
  if(NidLoop=="auto"){ # define the size of groups for each iterations of inla.run.many() calls
    # based on simulations, for simple to moderate models it is optimal to have a data size of ~12000
    # for complex models (~6+ likelihoods), it is optimal to have ~20000
    Nlik <- length(object$famLongi)+length(object$basRisk)
    # get an estimate of average data size per individual from fitted model:
    if(is.null(object$id)){
      ADS_i <- NsampleFE
    }else{
      ADS_i <- length(object$.args$data[[1]])*NsampleFE/length(na.omit(unique(object$.args$data[[object$id]])))
    }
    if(Nlik<6){
      NidLoop = round(12000 / ADS_i, 0)
    }else{
      NidLoop = round(20000 / ADS_i, 0)
    }
  }
  # baselineHaz = "smooth" | "interpolation"
  out <- NULL
  SumStats <- function(x) return(c(mean(x), sd(x), quantile(x, c(0.025, 0.5, 0.975))))
  if(!is.null(id)) idname <- id else idname <- object$id
  if(!is.null(object$id)) id <- object$id else if(is.null(id)) stop("Please specify individual id column name with argument 'id'")
  is_Long <- is_Surv <- FALSE
  idVect <- na.omit(unique(object$.args$data[[paste0("ID", object[["REstruc"]][[1]])]]))
  if(!as.character(object["REstruc"])=="NULL"){
    is_Long <- TRUE
  }
  if(!is.null(object$SurvInfo)){
    if(is.null(idVect) | length(idVect)==0){
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
        if(length(horizon)==1){
          if(horizon>max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]) & baselineHaz=="interpolation"){
            warning(paste0("The fitted model has baseline risk information up until value ",
                           max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]), " for survival outcome ", m, ". Since you ask for prediction at horizon ", horizon, " I will assume constant baseline hazard beyond the maximum available value. Alternatively, you can use baselineHaz='smooth' to use splines to predict the baseline hazard (for each sample). Alternatively, adding 'horizon' in the control options of the inla() call allows to extend the baseline beyond the last observed event time (linear extension based on last 2 values)."))
          }
        }else if(length(horizon)>1){
          if(T%in%c(horizon>max(object$.args$data[[paste0("baseline", m, ".hazard.values")]])) & baselineHaz=="interpolation"){
            warning(paste0("The fitted model has baseline risk information up until value ",
                           max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]), " for survival outcome ", m, ". Since you ask for prediction at horizon ", horizon[which(horizon>max(object$.args$data[[paste0("baseline", m, ".hazard.values")]]))[1]], " I will assume constant baseline hazard beyond the maximum available value. Alternatively, you can use baselineHaz='smooth' to use splines to predict the baseline hazard (for each sample). Alternatively, adding 'horizon' in the control options of the inla() call allows to extend the baseline beyond the last observed event time (linear extension based on last 2 values)."))
          }
        }
      }
    }
  }
  if(!is_Long & is_Surv){
    if(exists("newDataSurv") & is.null(newData)){
      newData <- newDataSurv
    }else if(exists("newData") & is.null(newDataSurv)){
      newDataSurv <- newData
    }
    if(is.null(object$timeVar) & !is.null(object$survOutcome)){
      object$timeVar <- as.character(object$SurvInfo[[1]]$nameTimeSurv)
    }
    if(!is.list(object$.args$data$Yjoint) & is.null(names(object$.args$data$Yjoint))){
      names(object$.args$data$Yjoint) <- "y1..coxph"
    }
  }
  if(inherits(newData, "tbl_df") || inherits(newData, "tbl")) {
    newData <- as.data.frame(newData)
  }
  if(!is_Long & !is_Surv) stop("Error, cannot recover ids from fitted model...")
  if(is_Surv & is.null(horizon)) stop("Please provide time horizon for prediction.")
  if(is_Surv & (length(horizon)>1 & length(horizon)!=length(unique(unique(newData[, object$id]))))) stop("Please provide either an unique horizon or a value for each id.")
  predL <- NULL
  predS <- NULL
  newPredS <- NULL
  if(is.null(object$id) & !is.null(id)) object$id <- id
  if(is.null(object$id)) stop("Please provide 'id' argument for new data.")
  ct <- object$misc$configs$contents
  if(is.null(ct)) stop("Please add argument 'cfg=TRUE' in control options when fitting the INLAjoint model to enable predictions.")
  if(ct$tag[1] == "Predictor") {
    ct$tag <- ct$tag[-1]
    ct$start <- ct$start[-1] - ct$start[2] + 1
    ct$length <- ct$length[-1]
  }
  paramVal <- object$misc$configs$config[[1]]$improved.mean # parameters value
  if(!is.null(startTime)) sTime <- startTime else sTime <- 0
  if(is.null(timePoints)){
    initTimePoints <- FALSE
    if(is.null(horizon)){
      timePoints <- seq(sTime, max(newData[, object$timeVar]), len=NtimePoints)
    }else if(length(horizon)==1){#} if(Csurv==0){
      timePoints <- seq(sTime, horizon, len=NtimePoints)
      # }else{
      #need to have a time point at Csurv there
    }
  }else{
    initTimePoints <- TRUE
  }
  firstID <- unique(newData[, object$id])[1]
  if(!silentMode) message("Sample...")
  SMPH <- INLA::inla.hyperpar.sample(NsampleHY, object)[rep(1:NsampleHY, each=NsampleFE),]
  SMP <- INLA::inla.rjmarginal(NsampleHY*NsampleFE, object)
  Nsample <- NsampleHY*NsampleFE

  if(FEonly){
    for(m in 1:M){
      idB_H <- c(sapply(1:M, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                          (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                             ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1)))
      SMP$samples[idB_H,] <- rowMeans(SMP$samples[idB_H,])
    }
  }
  if(!silentMode) message(paste0("Compute predictions..."))
  ParVal <- new("dgTMatrix", Dim=c(sum(ct$length), as.integer(Nsample)))
  ParValMode <- object$misc$configs$config[[1]]$improved.mean

  if(is_Surv){
    RWBH <- which(object$basRisk %in% c("rw1", "rw2"))
    if(length(RWBH)>0){
      SMP$samples[,1] <- ParValMode[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                                                 (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                                                    ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
                                      ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                       !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      ParVal[c(c(sapply(RWBH, function(x) ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]:
                          (ct$start[which(ct$tag %in% paste0("baseline", x, ".hazard"))]+
                             ct$length[which(ct$tag %in% paste0("baseline", x, ".hazard"))]-1))),
               ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                !ct$tag %in% paste0("baseline", 1:M, ".hazard"))]),] <- SMP$samples

    }else{
      SMP$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                                       !ct$tag %in% paste0("baseline", 1:M, ".hazard"))])]
      ParVal[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2) &
                                !ct$tag %in% paste0("baseline", 1:M, ".hazard"))]),] <- SMP$samples
    }
  }else{
    # use mode as first sample
    SMP$samples[,1] <- ParValMode[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))])]
    ParVal[c(ct$start[which(ct$tag %in% substr(rownames(SMP$samples), 1, nchar(rownames(SMP$samples))-2))]),] <- SMP$samples
  }
  nRE <- 0
  K <- 0 #number of longitudinal (written later, this is just to avoid errors when it is really 0)
  if(is_Long | !is.null(object[["REstrucS"]])){
    K <- length(object$famLongi) # number of longitudinal outcomes
    lenPV <- length(paramVal)
    SMPsel <- which(ct$length==1 &
                      substr(ct$tag, nchar(ct$tag)-2, nchar(ct$tag)-1)=="_L" |
                      substr(ct$tag, nchar(ct$tag)-3, nchar(ct$tag)-2)=="_L") # if >10 markers
    NamesH <- colnames(SMPH)
    nRE <- length(object[["REstruc"]])
    # Detect spatial effects and create filtered REstruc for SEL construction
    spatial_hyperpar_names <- unique(gsub("Phi for |Range for |Precision for |Stdev for ", "", colnames(SMPH)[sapply(colnames(SMPH), is_spatial_hyperpar)]))
    # Remove empty or invalid names
    spatial_hyperpar_names <- spatial_hyperpar_names[!is.na(spatial_hyperpar_names) & spatial_hyperpar_names != "" & !is.null(spatial_hyperpar_names)]
    # Create filtered REstruc and determine IID sampling strategy
    if(length(spatial_hyperpar_names) > 0 && !is.null(object[["REstruc"]]) && length(object[["REstruc"]]) > 0) {
      tryCatch({
        spatial_detection_results <- sapply(object[["REstruc"]], function(re) {
          if(is.null(re) || is.na(re) || re == "") return(FALSE)
          is_spatial_re(re, spatial_hyperpar_names)
        })
        if(!is.logical(spatial_detection_results)) {
          # no spatial effects detected
          spatial_re_indices <- integer(0)
        }else{
          spatial_re_indices <- which(spatial_detection_results)
        }
      }, error = function(e) {
        # If spatial detection fails, treat as no spatial effects
        spatial_re_indices <- integer(0)
      })
      if(length(spatial_re_indices) > 0) {
        # We have spatial effects - filter them out for INLA sampling
        filtered_re_struc <- object[["REstruc"]][-spatial_re_indices]
        nRE_for_SEL <- length(filtered_re_struc)
        all_effects_spatial <- (nRE_for_SEL == 0)
      }else{
        # No spatial effects found
        filtered_re_struc <- object[["REstruc"]]
        nRE_for_SEL <- nRE
        all_effects_spatial <- FALSE
      }
    }else{
      # No spatial hyperparameters detected
      filtered_re_struc <- object[["REstruc"]]
      nRE_for_SEL <- nRE
      all_effects_spatial <- FALSE
      spatial_re_indices <- integer(0)
    }
    if(is.null(object[["REstrucS"]])){
      if(nRE_for_SEL==1){
        BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE_for_SEL*Nsample) , as.integer(nRE_for_SEL*Nsample)))
        if(length(which(substr(colnames(SMPH), 1, 16)=="Precision for ID"))>0){
          if(length(which(substr(colnames(SMPH), 1, 16)=="Precision for ID"))>0){
            # Get precision indices excluding spatial effects
            precision_indices <- which(substr(colnames(SMPH), 1, 16)=="Precision for ID")
            filtered_indices <- precision_indices[!sapply(precision_indices, function(i) any(sapply(spatial_hyperpar_names, function(s) grepl(s, colnames(SMPH)[i]))))]
            if(length(filtered_indices) > 0) {
              diag(BD_Cmat) <- sqrt(1/SMPH[, filtered_indices])
            }
          }else if(length(which(substr(colnames(SMPH), 1, 12)=="Stdev for ID"))>0){ # SPDE
            # Filter out spatial SPDE effects
            stdev_indices <- which(substr(colnames(SMPH), 1, 12)=="Stdev for ID")
            filtered_stdev <- stdev_indices[!sapply(stdev_indices, function(i) any(sapply(spatial_hyperpar_names, function(s) grepl(s, colnames(SMPH)[i]))))]
            if(length(filtered_stdev) > 0) {
              diag(BD_Cmat) <- SMPH[, filtered_stdev]
            }
          }
        } # end nRE > 0 check
      }else if(nRE_for_SEL>1){
        # identify the position of the cholesky elements in hyperparameters
        if(object$corLong){
          PosH <- which(substr(NamesH, 1, 5)=="Theta" &
                          substr(NamesH, nchar(NamesH)-nchar(filtered_re_struc[1])-1,
                                 nchar(NamesH))==paste0("ID", filtered_re_struc[1]))
          # Block-Diagonal Cmatrix for all samples
          BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample)))
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
        }else{
          nRE_pk <- 1
          # Block-Diagonal Cmatrix for all samples (use filtered REs only)
          BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample), as.integer(nRE*Nsample)))
          for(k in 1:K){
            # Check if we have any filtered REs for this longitudinal component
            current_res_for_k <- filtered_re_struc[which(substr(filtered_re_struc, nchar(filtered_re_struc)-2, nchar(filtered_re_struc))==paste0("_L", k) |
                                                           substr(filtered_re_struc, nchar(filtered_re_struc)-3, nchar(filtered_re_struc))==paste0("_L", k))]

            if(length(current_res_for_k) == 0) {
              # No filtered REs for this component, skip
              next
            }

            if(nRE_pk <= length(filtered_re_struc)) {
              PosH <- which(substr(NamesH, 1, 5)=="Theta" &
                              substr(NamesH, nchar(NamesH)-nchar(filtered_re_struc[nRE_pk])-1,
                                     nchar(NamesH))==paste0("ID", filtered_re_struc[nRE_pk]))
            }else{
              PosH <- integer(0)  # No valid position found
            }

            nRE_k <- length(current_res_for_k)
            if(nRE_k==1){
              # Filter precision parameters to exclude spatial effects
              precision_matches <- which(substr(colnames(SMPH), 1, 16)=="Precision for ID" &
                                           (substr(colnames(SMPH), nchar(colnames(SMPH))-2, nchar(colnames(SMPH)))==paste0("_L", k) |
                                              substr(colnames(SMPH), nchar(colnames(SMPH))-3, nchar(colnames(SMPH)))==paste0("_L", k)))
              filtered_matches <- precision_matches[!sapply(precision_matches, function(i) any(sapply(spatial_hyperpar_names, function(s) grepl(s, colnames(SMPH)[i]))))]
              if(length(filtered_matches) > 0) {
                SMP_prec_k <- sqrt(1/SMPH[, filtered_matches])
              }else{
                # No matching precision parameters, create identity
                SMP_prec_k <- rep(1, Nsample)
              }
            }else{
              if(object$corRE[[1]]!=TRUE | length(object$corRE)>1){
                if(object$corRE[[k]]!=TRUE){
                  PosH <- sapply(filtered_re_struc, function(x) grep(x, NamesH))
                  if(length(PosH) > 0) {
                    SMP_prec_k <- SMPH[, PosH]
                  }else{
                    SMP_prec_k <- matrix(1, nrow=Nsample, ncol=nRE_k)
                  }
                }else{
                  if(length(PosH) > 0) {
                    L <- matrix(0, nrow=nRE_k, ncol=nRE_k)
                    # function to convert cholesky to precision
                    Chol_Prec <- function(x){
                      if(length(x[PosH]) >= nRE_k) {
                        diag(L) <- exp(x[PosH][1:nRE_k])
                        if(length(x[PosH]) > nRE_k) {
                          L[lower.tri(L)] <- x[PosH][-c(1:nRE_k)]
                        }
                        return(L %*% t(L))
                      }else{
                        return(diag(nRE_k))
                      }
                    }
                    SMP_prec_k <- apply(SMPH, 1, Chol_Prec)
                  }else{
                    SMP_prec_k <- array(diag(nRE_k), c(nRE_k, nRE_k, Nsample))
                  }
                }
              }else{
                if(length(PosH) > 0) {
                  L <- matrix(0, nrow=nRE_k, ncol=nRE_k)
                  # function to convert cholesky to precision
                  Chol_Prec <- function(x){
                    if(length(x[PosH]) >= nRE_k) {
                      diag(L) <- exp(x[PosH][1:nRE_k])
                      if(length(x[PosH]) > nRE_k) {
                        L[lower.tri(L)] <- x[PosH][-c(1:nRE_k)]
                      }
                      return(L %*% t(L))
                    }else{
                      return(diag(nRE_k))
                    }
                  }
                  SMP_prec_k <- apply(SMPH, 1, Chol_Prec)
                }else{
                  SMP_prec_k <- array(diag(nRE_k), c(nRE_k, nRE_k, Nsample))
                }
              }
            }
            # indices for BC_Cmat
            ind_BD_Cmat_k <- cbind(rep(rep(1:nRE_k, each=nRE_k), Nsample)+(rep(seq(nRE_pk, (Nsample*nRE), by=nRE), each=nRE_k^2)-1), rep(1:nRE_k, Nsample*nRE_k)+(rep(seq(nRE_pk, (Nsample*nRE), by=nRE), each=nRE_k^2)-1))
            if(object$corRE[[1]]!=TRUE | length(object$corRE)>1){
              if(object$corRE[[k]]!=TRUE){
                ind_BD_Cmat_k <- cbind(1:nrow(BD_Cmat), 1:nrow(BD_Cmat))
              }
            }
            # fill BD_Cmat
            BD_Cmat[ind_BD_Cmat_k] <- c(SMP_prec_k)
            nRE_pk <- nRE_pk + nRE_k # go to next block
          }
        }
      }
    }else{ # if there is at least a frailty, need to do the full model
      if(is_Long) nRES <- length(object[["REstrucS"]]) else nRES <- nRE
      if(nRE_for_SEL==1){ # only frailty
        BD_Cmat <- new("dgTMatrix", Dim=c(as.integer(nRE*Nsample) , as.integer(nRE*Nsample))) # adapt size
        # Find spatial effects and filter precision parameters
        spatial_effects <- unique(gsub("Phi for |Range for |Precision for |Stdev for ", "", colnames(SMPH)[sapply(colnames(SMPH), is_spatial_hyperpar)]))
        precision_matches <- which(substr(colnames(SMPH), 1, 16)=="Precision for ID")
        non_spatial_matches <- precision_matches[!sapply(precision_matches, function(i) any(sapply(spatial_effects, function(s) grepl(s, colnames(SMPH)[i]))))]
        if(length(non_spatial_matches) > 0) {
          diag(BD_Cmat) <- sqrt(1/SMPH[, non_spatial_matches])
        }
      }else if(nRE_for_SEL>1){
        # need to do the full model to get posteriors to sample from as there is at least longi and frailty here
      }
    }
    ResErrFixed <- vector("list", K)
    if(is.null(object[["REstrucS"]])){
      # Filter REnames to exclude spatial effects
      REnames <- c(sapply(filtered_re_struc, function(x) paste0("ID", x)))
    }else{
      if(!as.character(object["REstruc"])=="NULL"){
        # Filter REnames to exclude spatial effects
        REnames <- c(sapply(filtered_re_struc, function(x) paste0("ID", x)))
        REnamesS <- object[["REstrucS"]]
      }else{
        REnames <- REnamesS <- object[["REstrucS"]]
      }
    }
    posRE <- sort(ct$start[unlist(sapply(REnames, function(x) which(ct$tag==x)))])
    ordRE <- order(order(ct$start[unlist(sapply(REnames, function(x) which(ct$tag==x)))]))
    assocNs <- object$assoc
    assocNa <- object$assoc_Names
    if(is_Surv){
      assocPos <- unlist(sapply(assocNs, function(x) grep(x, ct$tag)))
      # identify the longitudinal needed for association
      # first identify shared part from longitudinal (no duplicates, so if CV from longitudinal 1 is shared twice, we need to repeat it)
      OutcNam <- substr(names(object$.args$data$Yjoint), 1, nchar(names(object$.args$data$Yjoint))-1)
      if(!is.null(names(object$.args$data$Yjoint))){
        outcomeAssoc <- names(object$.args$data$Yjoint)[unlist(sapply(1:length(OutcNam), function(x) if(length(grep(OutcNam[x], ct$tag))!=0) return(x)))]
        outcomeAssoc2 <- sapply(outcomeAssoc, function(x) strsplit(x, split = "_S")[[1]][1])
        requiredAssoc <- sapply(assocNs, function(x) strsplit(x, split = "_S")[[1]][1])
        patternAsso <- unname(sapply(requiredAssoc, function(x) which(x==outcomeAssoc2)))
      }else{
        outcomeAssoc <- NULL
        outcomeAssoc2 <- NULL
        requiredAssoc <- NULL
        patternAsso <- NULL
      }
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
        SRE_inda <- c(unlist(sapply(paste0("SRE_", object[["REstruc"]]), function(x) grep(x,assocNs))))
        SRE_inda2 <- c(unlist(sapply(assocNs, function(x) which(sapply(substr(REnames, 3, nchar(REnames)), function(xx) grep(xx, x))==1))))
      }else{
        SRE_inda <- NULL
        SRE_inda2 <- NULL
      }
      if(length(SRE_inda)>0){
        patternAsso2 <- unname(unlist(patternAsso[-SRE_inda]))
        patternAsso <- 1:length(patternAsso)
        assocNs2 <- assocNs[-SRE_inda]
      }else{
        patternAsso2 <- patternAsso
        assocNs2 <- assocNs
      }
    }
  }
  if(idLoop) idLoopSet <- FALSE else idLoopSet <- TRUE # used when all individuals random effects done in 1 call
  RErun_iter <- 0
  newRErun <- NULL
  reloadCT <- TRUE
  horizonF <- horizon # keep it when horizon is a vector
  for(idPred in unique(newData[, object$id])){
    ct2 <- ct
    if(length(horizonF)>1){ # horizon is different for each id
      horizon <- horizonF[which(unique(newData[, object$id])==idPred)]
      if(initTimePoints) NtimePoints <- length(timePoints)
      if(!initTimePoints){
        timePoints <- seq(sTime, horizon, len=NtimePoints)
      }
    }
    if(NidLoop=="auto"){
      NidLoop <- 1
    }
    # split data to loop over groups of individuals
    ND_split <- split(unique(newData[, object$id]), ceiling(seq_along(unique(newData[, object$id]))/NidLoop))
    ND_id <- unname(which(sapply(ND_split, function(x) idPred %in% x)))
    curID <- ND_split[[ND_id]]
    if(RErun_iter < ND_id){# need to estimate new RE posteriors?
      newRErun <- TRUE
      RECOUNT_ <- 1
    }else{
      newRErun <- FALSE
      RECOUNT_ <- RECOUNT_ + 1
    }
    RErun_iter <- ND_id
    ND <- newData[newData[, object$id,] %in% ND_split[[ND_id]],,drop=FALSE]
    idPredt <- which(unique(ND[, object$id])==idPred)
    ND[, object$id] <- sapply(ND[, object$id], function(x) (1:length(unique(ND[, object$id])))[which(unique(ND[, object$id])==x)])
    if(!is.null(object$lonFacChar) & length(which(names(object$lonFacChar) %in% colnames(ND)))>0){
      for(Fi in which(names(object$lonFacChar) %in% colnames(ND))){
        # ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])] <- factor(gsub(" ","", gsub("[^[:alnum:] ]","", ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])])), levels=gsub(" ","", gsub("[^[:alnum:] ]","", object$lonFacChar[[Fi]])))
        ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])] <- factor(gsub(" ","", gsub("[^[:alnum:] ]","", ND[, which(colnames(ND)==names(object$lonFacChar)[Fi])])), levels=object$lonFacChar[[Fi]])
      }
    }
    if(is_Long & is.null(Csurv)){
      if(object$timeVar %in% colnames(ND)){
        TPO <- sort(unique(c(timePoints, max(ND[, object$timeVar]))))
      }else{
        TPO <- timePoints
        ND <- cbind(ND, 0)
        colnames(ND)[length(colnames(ND))] <- object$timeVar
      }
      NTP <- length(TPO)
    }else if(!is_Long & is.null(Csurv)){
      TPO <- timePoints
      NTP <- NtimePoints
    }else if(Csurv==0){
      TPO <- timePoints
      NTP <- NtimePoints
    }else{
      TPO <- sort(unique(c(timePoints, Csurv)))
      NTP <- length(TPO)#NtimePoints+1
    }
    call.new2 <- object$call
    TXT1 <- NULL
    if(is_Surv){
      if(is_Long & !is.null(newDataSurv)){
        if(NidLoop!=FALSE){
          NDS <- newDataSurv[newDataSurv[, object$id] %in% ND_split[[ND_id]],]
          NDS[, object$id] <- sapply(NDS[, object$id], function(x) (1:length(unique(NDS[, object$id])))[which(unique(NDS[, object$id])==x)])
        }else{
          NDS <- newDataSurv
        }
      }else{
        if(!is.null(object[["REstrucS"]])){
          NDS <- ND
        }else{
          NDS <- ND[which(!duplicated(ND[,id], fromLast=TRUE)),, drop=FALSE]
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
          if(!S_Outc %in% colnames(NDS)){ # condition on specific time
            if(!is.null(Csurv)){
              NDS <- cbind(NDS, 0)
              colnames(NDS)[length(colnames(NDS))] <- S_Outc
            }else if(!is.null(newDataSurv) & as.character(object$SurvInfo[[1]]$survOutcome) %in% colnames(newDataSurv)){ # condition on newDataSurv
              NDS <- cbind(NDS, newDataSurv[, as.character(object$SurvInfo[[1]]$survOutcome)])
              colnames(NDS)[length(colnames(NDS))] <- S_Outc
            }else{
              NDS <- cbind(NDS, 0)
              colnames(NDS)[length(colnames(NDS))] <- S_Outc
            }
          }
          if(!S_nam %in% colnames(NDS)){
            if(is.null(object$timeVar)){
              colTS <- which(colnames(ND) %in% unlist(sapply(object$SurvInfo, function(x) x$nameTimeSurv)))
              mTS <- max(ND[colTS])
            }else{
              if(!is.null(newDataSurv) & as.character(object$SurvInfo[[1]]$nameTimeSurv) %in% colnames(newDataSurv)){
                mTS <- newDataSurv[, as.character(object$SurvInfo[[1]]$nameTimeSurv)]
              }else if(object$timeVar %in% colnames(ND)){
                if(!is.null(Csurv)){
                  mTS <- Csurv
                }else{
                  mTS <- ND[!duplicated(ND[[object$id]], fromLast = T), object$timeVar]
                }
              }else{
                mTS <- 0
              }
            }
            NDS <- cbind(NDS, mTS)
            colnames(NDS)[length(colnames(NDS))] <- S_nam
            if(TRUE %in% (mTS!=0)){
              ND <- cbind(ND, mTS[ND[[object$id]]])
            }else{
              ND <- cbind(ND, 0)
            }
            colnames(ND)[length(colnames(ND))] <- S_nam
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
      if(!is.null(object$lonFacChar) & length(which(names(object$lonFacChar) %in% colnames(NDS)))>0){
        for(Fi in which(names(object$lonFacChar) %in% colnames(NDS))){
          # NDS[, which(colnames(NDS)==names(object$lonFacChar)[Fi])] <- factor(gsub(" ","", gsub("[^[:alnum:] ]","", NDS[, which(colnames(NDS)==names(object$lonFacChar)[Fi])])), levels=gsub(" ","", gsub("[^[:alnum:] ]","", object$lonFacChar[[Fi]])))
          NDS[, which(colnames(NDS)==names(object$lonFacChar)[Fi])] <- factor(gsub(" ","", gsub("[^[:alnum:] ]","", NDS[, which(colnames(NDS)==names(object$lonFacChar)[Fi])])), levels=object$lonFacChar[[Fi]])
        }
      }
      if(max(ND[object$timeVar])>horizon & is.null(Csurv)) warning(paste0("horizon = ", horizon, " and there are observations up to ", max(ND[object$timeVar]), ". It is likely not what you want but you can use Csurv argument if you want to to force predictions conditional on future observations."))
      horizon2 <- max(TPO)+0.0001
      # SdataPred <- ND[!duplicated(ND[, object$id]),]
      if(!is.null(object[["REstrucS"]])){
        SdataPred <- NDS[NDS[,id]==idPredt,]
      }else{
        SdataPred <- NDS[NDS[,id]==idPredt,][1,]
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
          # if(SdataPred[nrow(SdataPred), SVO]==1 & !is.null(object[["REstrucS"]])){
          #   SdataPred <- SdataPred[c(1:nrow(SdataPred), nrow(SdataPred)),]
          # }
          SdataPred[nrow(SdataPred), SVO] <- 0
        }
        # remove truncation to force predict at all given time points
        if(!is.null(object$SurvInfo[[m]]$nameTrunc) & !is.null(startTime)){
          SdataPred[, which(colnames(SdataPred)==object$SurvInfo[[m]]$nameTrunc)] <- min(TPO)
        }else{
          startTime <- 0
        }
      }
      if(!is.null(object$dataSurv)){
        if(paste0(object$dataSurv)[1]=="list"){
          if(length(object[["REstrucS"]])>9) stop("Predictions not implemented for 10+ frailties, contact INLAjoint@gmail.com")
          for(m in 1:(length(paste0(object$dataSurv))-1)){ # all lines (changed to only last line as this is design (uData has all lines))
            if(length(grep(paste0("_S", m), substr(object[["REstrucS"]],
                                                   start=nchar(object[["REstrucS"]])-2,
                                                   stop=nchar(object[["REstrucS"]]))))>0){
              assign(paste0(object$dataSurv)[m+1], SdataPred[nrow(SdataPred),])
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
    }else{
      NDS <- NULL
      M <- 0
    }
    if(is_Long){
      LdataPred <- ND[ND[,id]==idPredt,][rep(1, length(TPO)), ]
      LdataPred[, object$timeVar] <- TPO
      if(is.null(startTime)){ # start predictions from first observed time
        TPO <- TPO[which(TPO>=min(ND[which(ND[, object$id]==idPredt), object$timeVar]))]
        NTP <- length(TPO)
        LdataPred <- LdataPred[which(LdataPred[, object$timeVar]>=min(ND[which(ND[, object$id]==idPredt), object$timeVar])),]
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
    Yjoint_original <- NEWdata$Yjoint
    survPart <- NULL
    # if(is_Surv & M>1 & K==0) survPart <- unique(c(sapply(NEWdata$Yjoint, function(x) which(!is.na(x))))) # this could be simplified, all points are included in this case
    # if(is_Surv & (M+K)>1 & K>0) survPart <- c(unlist(sapply(1:M, function(x) which(!is.na(eval(parse(text=paste0("NEWdata$Yjoint$y", x, "..coxph"))))))))
    # if(is_Surv & !is_Long & (M==1)) survPart <- c(unlist(which(!is.na(eval(parse(text=paste0("NEWdata$Yjoint")))))))
    if(is_Surv){
      for(m in 1:M){
        survPart <- c(survPart, which(!is.na(eval(parse(text=paste0("NEWdata$y", m, "..coxph"))))))
      }
    }
    if(!is.list(NEWdata)) NEWdata <- as.list(as.data.frame(NEWdata))
    if(is_Long | !is.null(object[["REstrucS"]])){
      ###              ###
      ### LONGITUDINAL ###
      ###              ###
      ctL <- ct
      inND <- as.integer(ND[, object$id])
      ND <- ND[rep(1:nrow(ND), NsampleFE),]
      ND[, object$id] <- rep(inND, NsampleFE) + rep(max(inND), NsampleFE*length(inND))*rep(0:(NsampleFE-1), each=length(inND))
      if(is_Surv){
        inNDS <- as.integer(NDS[, object$id])
        NDS <- NDS[rep(1:nrow(NDS), NsampleFE),]
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
      if(!is.null(object[["REstrucS"]]) | is_Surv){
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
            assign(paste0(object$dataSurv), NDS)
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

        # Filter out spatial effects from formula if set.samples is provided
        if(!is.null(set.samples)){
          # Split formula into terms
          formula_terms <- strsplit(FRM2, " \\+ ")[[1]]
          # Keep only non-spatial terms
          filtered_terms <- formula_terms[!sapply(formula_terms, function(term) {
            any(sapply(names(set.samples), function(spatial_name) grepl(spatial_name, term)))
          })]
          FRM2 <- paste(filtered_terms, collapse = " + ")
        }
        SPLIT_n <- strsplit(FRM2, " n = (.*?),")[[1]] # change the length of iid random effects as data is different
        # Filter out spatial effects from SPLIT_n if set.samples is provided
        if(!is.null(set.samples)) {
          # Keep only non-spatial terms in SPLIT_n
          spatial_indices <- sapply(SPLIT_n, function(term) {
            any(sapply(names(set.samples), function(spatial_name) grepl(spatial_name, term)))
          })
          filtered_split <- SPLIT_n[!spatial_indices]
          # Also filter nre_prT to match the filtered SPLIT_n
          if(exists("nre_prT") && length(nre_prT) > 0) {
            # Remove nre_prT entries corresponding to filtered SPLIT_n entries
            # Note: nre_prT corresponds to SPLIT_n[1:(length(nre_prT))]
            if(length(nre_prT) <= length(SPLIT_n)) {
              spatial_nre_indices <- spatial_indices[1:length(nre_prT)]
              nre_prT <- nre_prT[!spatial_nre_indices]
            }
          }
          SPLIT_n <- filtered_split
        }

        # recover length of each iid random effect groups
        nre_pr <- NULL
        if(is_Long & is.null(object[["REstrucS"]])){
          if(object$corLong) rmvCL = 0 else rmvCL=1
          if(object$corRE[[1]]) if(length(object$famLongi) != (length(SPLIT_n)-rmvCL) & rmvCL==1) if(length(object$famLongi) != length(SPLIT_n)) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
          if(length(SPLIT_n) !=2  & rmvCL==0) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
          for(nre_p in 1:length(object$famLongi)){
            if(nre_p<10){
              nre_10p = 0
            }else if(nre_p>=10){
              nre_10p = 1
            }
            if(!is.null(object$corRE)){
              if(object$corRE[[1]]!=TRUE | length(object$corRE)>1){
                if(object$corRE[[nre_p]]!=TRUE){
                  nre_pr <- c(nre_pr, length(unique(ND[,id])))
                }else{
                  nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
                }
              }else{
                nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
              }
            }else{
              nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
            }
          }
          if(object$corLong) nre_prT <- sum(nre_pr) else nre_prT <- nre_pr

          # Adjust nre_prT for spatial filtering: when spatial effects are removed,
          # the IID random effect should use the number of unique IDs in prediction data
          if(!is.null(set.samples)) {
            n_unique_ids <- length(unique(ND[,id]))
            if(object$corLong) nre_prT <- n_unique_ids else nre_prT <- rep(n_unique_ids, length(nre_pr))
          }
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
            nre_pr <- c(nre_pr, length(grep(paste0("_L", nre_p), substr(object[["REstruc"]], start=nchar(object[["REstruc"]])-2-nre_10p, stop=nchar(object[["REstruc"]]))))*length(unique(ND[,id])))
          }
          if(object$corLong) nre_prT <- sum(nre_pr) else nre_prT <- nre_pr
          if(!is.null(set.samples)) {
            n_unique_ids <- length(unique(ND[,id]))
            if(object$corLong) nre_prT <- n_unique_ids else nre_prT <- rep(n_unique_ids, length(nre_pr))
          }
          for(nre_p in 1:length(object[["REstrucS"]])){ # then surv frailty random effects
            if(nre_p<10){
              nre_10p = 0
            }else if(nre_p>=10){
              nre_10p = 1
            }
            nre_pr <- c(length(grep(paste0("_S", nre_p), substr(object[["REstrucS"]], start=nchar(object[["REstrucS"]])-2-nre_10p, stop=nchar(object[["REstrucS"]]))))*length(unique(ND[,id])), nre_pr)
          }
        }else if(!is_Long & !is.null(object[["REstrucS"]])){
          if(length(object[["REstrucS"]]) != (length(SPLIT_n)-1) & length(SPLIT_n)>1) stop("I found a mismatch for some internal computations, please report to INLAjoint@gmail.com")
          for(nre_p in 1:length(object[["REstrucS"]])){
            if(nre_p<10){
              nre_10p = 0
            }else if(nre_p>=10){
              nre_10p = 1
            }
            nre_pr <- c(nre_pr, length(grep(paste0("_S", nre_p), substr(object[["REstrucS"]], start=nchar(object[["REstrucS"]])-2-nre_10p, stop=nchar(object[["REstrucS"]]))))*length(unique(ND[,id])))
          }
          nre_prT <- nre_pr
        }
        if(object$corLong){
          FRM3 <- paste(paste(sapply(1:(1+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
        }else{
          if((length(object$famLongi)+length(object[["REstrucS"]]))>0 & length(SPLIT_n)>1){
            if(!is.null(object$corRE)){
              if(object$corRE[[1]]!=TRUE | length(object$corRE)>1){
                for(k in 1:length(object$corRE)){
                  if(object$corRE[[k]]!=TRUE){
                    FRM3 <- paste(paste(sapply(1:nRE_k, function(x) paste0(SPLIT_n[x], " n = ", nre_prT[1], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
                  }else{
                    FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
                  }
                }
              }else{
                FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
              }
            }else{
              FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
            }
          }else{
            FRM3 <- FRM2
          }
        }
        # remove baseline from formula as it comes from full model (sampled as fixed effect)
        if(is_Surv){
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

        # Filter out spatial effects from formula if set.samples is provided
        if(!is.null(set.samples)) {
          # Split formula into terms
          formula_terms <- strsplit(FRM2, " \\+ ")[[1]]
          # Keep only non-spatial terms
          filtered_terms <- formula_terms[!sapply(formula_terms, function(term) {
            any(sapply(names(set.samples), function(spatial_name) grepl(spatial_name, term)))
          })]
          FRM2 <- paste(filtered_terms, collapse = " + ")
        }

        SPLIT_n <- strsplit(FRM2, " n = (.*?),")[[1]] # change the length of iid random effects as data is different

        # Filter out spatial effects from SPLIT_n if set.samples is provided
        if(!is.null(set.samples)) {
          # Keep only non-spatial terms in SPLIT_n
          spatial_indices <- sapply(SPLIT_n, function(term) {
            any(sapply(names(set.samples), function(spatial_name) grepl(spatial_name, term)))
          })
          filtered_split <- SPLIT_n[!spatial_indices]

          # Also filter nre_prT to match the filtered SPLIT_n
          if(exists("nre_prT") && length(nre_prT) > 0) {
            # Remove nre_prT entries corresponding to filtered SPLIT_n entries
            # Note: nre_prT corresponds to SPLIT_n[1:(length(nre_prT))]
            if(length(nre_prT) <= length(SPLIT_n)) {
              spatial_nre_indices <- spatial_indices[1:length(nre_prT)]
              nre_prT <- nre_prT[!spatial_nre_indices]
            }
          }

          SPLIT_n <- filtered_split
        }

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
          if(object$corLong) nre_prT <- sum(nre_pr) else nre_prT <- nre_pr

          # Apply spatial filtering to nre_prT in the second path as well
          if(!is.null(set.samples)) {
            # Adjust nre_prT for spatial filtering: use number of unique IDs in prediction data
            n_unique_ids <- length(unique(ND[,id]))
            if(object$corLong) nre_prT <- n_unique_ids else nre_prT <- rep(n_unique_ids, length(nre_pr))
          }
        }
        if(object$corLong){
          FRM3 <- paste(paste(sapply(1:(1+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
        }else{
          FRM3 <- paste(paste(sapply(1:(length(object$famLongi)+length(object[["REstrucS"]])), function(x) paste0(SPLIT_n[x], " n = ", nre_prT[x], ","), simplify=F), collapse=''), SPLIT_n[length(SPLIT_n)], collapse='')
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
      if(!is.list(uData)) uData <- as.list(as.data.frame(uData))
      nL_K <- length(uData[[1]])
      # now we prepare the precision matrix for all samples (large block diagonal matrix)
      IDshift <- 0
      # A matrix for offset computation
      # ids to select the elements to keep in latent part of samples
      # baseline => substr(ct$tag, 1, 8)=="baseline" |
      A_off <- new("dgTMatrix", Dim=c(nL_K, sum(ct$length)))
      if(is_Long) A_off[, ct$start[SMPsel]] <- do.call(cbind, sapply(uData[ct$tag[SMPsel]], function(x) replace(x, is.na(x), 0), simplify=F))
      if(!is.null(object[["REstrucS"]]) | is_Surv){
        SMPselS <- which(ct$length==1 &
                           substr(ct$tag, nchar(ct$tag)-2, nchar(ct$tag)-1)=="_S" |
                           substr(ct$tag, nchar(ct$tag)-3, nchar(ct$tag)-2)=="_S") # if >10 markers
        if(length(SMPselS)>0){
          A_off[, ct$start[SMPselS]] <- do.call(cbind, sapply(uData[ct$tag[SMPselS]], function(x) replace(x, is.na(x), 0), simplify=F))
        }
      }
      if(!is.null(set.samples)){
        for(rsmp in 1:length(set.samples)){
          Nrsmp <- names(set.samples)[rsmp]
          id_L <- as.integer(strsplit(Nrsmp, "_L")[[1]][2])
          id_S <- as.integer(strsplit(Nrsmp, "_S")[[1]][2])
          if(!is.na(id_L)){
            colSEL <- 1:nre_prT + nre_prT*(id_L-1)
          }else if(!is.na(id_S)){
            colSEL <- which(!is.na(uData[[paste0("baseline", id_S, ".hazard.idx")]]))
          }
          A_off[colSEL, ct$start[ct$tag==Nrsmp]] <- 1
          if(is.null(dim(set.samples[[rsmp]]))){
            ParVal[ct$start[ct$tag==Nrsmp], ] <- set.samples[[rsmp]]
            ParValMode[ct$start[ct$tag==Nrsmp]] <- mean(set.samples[[rsmp]])
          }else{
            ParVal[ct$start[ct$tag==Nrsmp], ] <- set.samples[[rsmp]][idPredt, ]
            ParValMode[ct$start[ct$tag==Nrsmp]] <- mean(set.samples[[rsmp]][idPredt, ])
          }
        }
        # need to remove the corresponding random effect from formula
        # if no RE left => skip inla call
        if(length(c(object[["REstruc"]], object[["REstrucS"]]))==1){
          if(names(set.samples) == c(object[["REstruc"]], object[["REstrucS"]])){
            newRErun <- FALSE # skip inla() call as the unique RE is pre-sampled
          }
        }
        # add shared and scaled random effect
        if(exists("Nrsmp")){ # if some samples are set, they may need scaling for shared frailty
          if(length(unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))>0){
            for(ias in 1:length(unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))){
              m_inti <- unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH)))))[ias]
              m_intiCT <- which(ct2$tag == gsub("Beta for ", "", colnames(SMPH)[m_inti]))
              id_S <- strsplit(ct2$tag[m_intiCT], paste0(gsub("^ID", "", Nrsmp), "_S"))[[1]][2]
              A_off[which(!is.na(uData[[paste0("baseline", id_S, ".hazard.idx")]])), ct2$start[m_intiCT]] <- 1
              # compute scaled frailty term and insert in shared part
              if(is.null(dim(set.samples[[ias]]))){
                ParVal[ct2$start[m_intiCT], ] <- set.samples[[ias]] * SMPH[, m_inti]
                ParValMode[ct2$start[m_intiCT]] <- mean(set.samples[[ias]] * SMPH[, m_inti])
              }else{
                ParVal[ct2$start[m_intiCT], ] <- set.samples[[ias]][idPredt,] * SMPH[, m_inti]
                ParValMode[ct2$start[m_intiCT]] <- mean(set.samples[[ias]][idPredt,] * SMPH[, m_inti])
              }
            }
          }
        }
      }
      # set baseline in A_off
      if(is_Surv){
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
      if(!is.null(assocNa)){
        if(Nsample==1){
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
          offSet <- A_off %*% ParVal[, 1:(Nsample)] # samples
          for(a_id in 1:length(assocNa)){
            # grab values of linear combination of fixed effects to share
            LP_sh <- offSet[which(!is.na(uData$Yjoint[[grep(assocNa[a_id], names(uData$Yjoint))]])),]
            # scale it
            # LP_shsc <- LP_sh * c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
            #                      SMPH[1:(Nsample-1), grep(assocNs[a_id], colnames(SMPH))])
            if(!is.null(dim(LP_sh))){
              LP_shsc <- sapply(1:length(c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                           SMPH[1:(Nsample-1), grep(assocNs[a_id], colnames(SMPH))])),
                                function(x) LP_sh[, x]*c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                                         SMPH[1:(Nsample-1), grep(assocNs[a_id], colnames(SMPH))])[x])
            }else{
              LP_shsc <- sapply(1:length(c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                           SMPH[1:(Nsample-1), grep(assocNs[a_id], colnames(SMPH))])),
                                function(x) LP_sh[x]*c(object$summary.hyperpar[grep(assocNs[a_id], rownames(object$summary.hyperpar)), "0.5quant"],
                                                       SMPH[1:(Nsample-1), grep(assocNs[a_id], colnames(SMPH))])[x])
            }
            # add it to offset
            LPS_index <- which(!is.na(uData[[grep(paste0("^", assocNa[a_id], "$"), names(uData))]]))
            offSet[LPS_index,] <- offSet[LPS_index, ] + LP_shsc
          }
        }
      }else{
        offSet <- A_off %*% ParVal[, 1:(Nsample)]
      }
      if(is_Long){
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
      if(newRErun){ # full inla call
        RMNk <- object$.args$control.fixed$remove.names
        object$.args$control.fixed$remove.names <- c(object$.args$control.fixed$remove.names, rownames(object$summary.fixed))
        SEL <- NULL
        re_SUM <- 0
        if(!is.null(nre_pr)){
          for(re_i in 1:length(nre_pr)){
            if(!is.null(object[["REstrucS"]])){
              if(re_i<=length(object[["REstrucS"]])){
                SEL <- append(SEL, list((1:length(unique(ND[,id])))))
              }else{
                for(re_j in 1:length(filtered_re_struc)){
                  SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_j-1)))
                }
              }
            }else{
              # Only process if there are matching effects in filtered_re_struc
              matching_effects <- grep(paste0("_L", re_i), filtered_re_struc)
              if(length(matching_effects) > 0) {
                for(re_j in 1:length(matching_effects)){
                  # for(re_j in 1:nre_pr[re_i]){
                  if(!object$corLong){
                    if(object$corRE[[1]]!=TRUE | length(object$corRE)>1){
                      if(object$corRE[[k]]!=TRUE){
                        SEL <- append(SEL, list((1:length(unique(ND[,id])))))
                      }else{
                        SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_j-1)))
                      }
                    }else{
                      SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_j-1)))
                    }
                  }else{
                    if(re_j>1) re_SUM <- re_SUM+1
                    SEL <- append(SEL, list((1:length(unique(ND[,id])))+length(unique(ND[,id]))*(re_i+re_SUM-1)))
                    # SEL <- append(SEL, list(1:nre_pr))
                  }
                }
              }
            }
          }
        }
        if(!is.null(object[["REstruc"]])){
          if(length(spatial_re_indices) > 0){
            names_reL <-paste0("ID", object[["REstruc"]][-spatial_re_indices])
          }else{
            names_reL <-paste0("ID", object[["REstruc"]])
          }
        }else{
          names_reL <- NULL
        }
        # if(object$corLong) names_reL <- names_reL[1]
        if(!is.null(object[["REstrucS"]])) names_reS <- object[["REstrucS"]] else names_reS <- NULL
        names(SEL) <- c(names_reS, names_reL)
        # SPATIAL FILTERING: Filter SEL to exclude spatial effects from posterior sampling
        if(length(spatial_hyperpar_names) > 0 && length(names(SEL)) > 0) {
          # Identify spatial effects in SEL names
          spatial_sel_indices <- which(sapply(names(SEL), function(name) {
            # Check if this name corresponds to a spatial effect
            any(sapply(spatial_hyperpar_names, function(spatial_name) {
              grepl(spatial_name, name, fixed=TRUE)
            }))
          }))
          if(length(spatial_sel_indices) > 0) {
            SEL <- SEL[-spatial_sel_indices]
            # Also filter names_reS and names_reL to match the filtered SEL
            if(!is.null(names_reS)) {
              spatial_reS_indices <- which(sapply(names_reS, function(name) {
                any(sapply(spatial_hyperpar_names, function(spatial_name) {
                  grepl(spatial_name, name, fixed=TRUE)
                }))
              }))
              if(length(spatial_reS_indices) > 0) {
                names_reS <- names_reS[-spatial_reS_indices]
              }
            }
            if(!is.null(names_reL)) {
              spatial_reL_indices <- which(sapply(names_reL, function(name) {
                # Remove "ID" prefix for detection
                clean_name <- gsub("^ID", "", name)
                any(sapply(spatial_hyperpar_names, function(spatial_name) {
                  grepl(spatial_name, clean_name, fixed=TRUE)
                }))
              }))
              if(length(spatial_reL_indices) > 0) {
                names_reL <- names_reL[-spatial_reL_indices]
              }
            }
          }
        }
        infoBHSEL <- 0
        if(Nsample>1) RE_values <- NULL
        if(is_Surv) BH_values <- NULL
        if(Nsample==1){
          if(is_Surv){
            TETA <- data.frame(object$misc$theta.mode[-grep("baseline", object$misc$theta.tags)])
          }else{
            TETA <- data.frame(object$misc$theta.mode)
          }
          # OFFSET <- data.frame(offS)
          uData <- append(uData, list("off"=as.matrix(offS)))
        }else{
          if(TRUE %in% sapply(c("rw1", "rw2"), function(x) x %in% unlist(object$basRisk)) & is_Surv){
            TETA <- sapply(1:Nsample, function(S) sapply(1:length(unname(SMPH[S,])), function(x) object$misc$to.theta[[x]](unname(SMPH[S,,drop=FALSE])[x]))[-grep("baseline", object$misc$theta.tags)])
            if(is.null(dim(TETA))) TETA <- t(data.frame(TETA))
          }else{
            TETA <- sapply(1:Nsample, function(S) sapply(1:length(unname(SMPH[S,])), function(x) object$misc$to.theta[[x]](unname(SMPH[S,,drop=FALSE])[x])))
            if(is.null(dim(TETA))) TETA <- t(data.frame(TETA))
          }
          uData <- append(uData, list("off"=as.matrix(offSet)))
        }
        TETA <- TETA[, seq(1, NsampleHY*NsampleFE, NsampleFE),drop=FALSE]
        offS_NEW <- matrix(NA, nrow = length(uData[[1]]), ncol=NsampleHY)
        for(n_HY in 1:NsampleHY){
          offS_HY <- uData$off[ ,1:NsampleFE + rep((n_HY-1)*NsampleFE, NsampleFE)]
          if(is_Surv){
            for(m in 1:M){
              # set survival part of offset for each Hyperpar sample (NsampleHY)
              # and each FE sample (NsampleFE)
              BH_m <- which(!is.na(uData[[paste0("expand", m, "..coxph")]]))
              n_FEi <- length(BH_m)/NsampleFE
              offS_NEW[BH_m, n_HY] <- c(sapply(1:NsampleFE, function(x) offS_HY[BH_m[1:n_FEi], x]))
            }
          }
          # LP_K <- sapply(object$longOutcome, function(x) grep(x, names(uData$Yjoint)))
          LP_K <- sapply(1:length(object$longOutcome), function(x) which(names(uData$Yjoint)==paste0(object$longOutcome[[x]], "_L", x)))
          if(is_Long){
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
        # INLA::inla.tempdir()
        INLA::inla.setOption(malloc.lib='compiler')
        INLA::inla.setOption(INLAjoint.features=TRUE)
        object$.args$control.inla$compute.initial.values=FALSE
        wd <- .inla_tempdir_safe()#"model.files"
        # unlink(wd, recursive = TRUE)
        if(length(which(object$.args$control.predictor$link!=1))>0) warning("Link function is not default, this has to be added here and has not yet been done. Please contact INLAjoint@gmail.com")
        if(!silentMode & REmsg) message("Estimate conditional posterior of random effects (N = ", length(unique(newData[, object$id])), ")...")
        if(REmsg) REmsg <- FALSE
        if(length(unique(newData[, object$id]))>=length(curID) & !silentMode & NidLoop>1) message(paste0("... id ", curID[1], " to ", tail(curID, 1), "..."))
        if(length(unique(newData[, object$id]))>=length(curID) & !silentMode & NidLoop==1) message(paste0("... id ", curID[1], "..."))
        if(Nsample==1) TETA <- TETA[,1,drop=FALSE]
        # Filter TETA to remove spatial hyperparameters when set.samples is provided
        if(!is.null(set.samples)) {
          # Find indices of spatial hyperparameters in theta.tags
          spatial_theta_indices <- c()
          for(spatial_name in names(set.samples)) {
            # Look for hyperparameters containing the spatial effect name
            # Updated patterns based on actual theta.tags format
            spatial_patterns <- c(paste0("Log precision for ", spatial_name),
                                  paste0("Logit phi for ", spatial_name),
                                  paste0("Phi for ", spatial_name),
                                  paste0("Range for ", spatial_name),
                                  paste0("Precision for ", spatial_name),
                                  paste0("Beta_intern for SRE_", gsub("^ID", "", spatial_name)))
            for(pattern in spatial_patterns) {
              matches <- grep(pattern, object$misc$theta.tags, fixed = TRUE)
              spatial_theta_indices <- c(spatial_theta_indices, matches)
            }
            # Also search for any theta.tags containing the spatial name
            matches <- grep(spatial_name, object$misc$theta.tags, fixed = TRUE)
            spatial_theta_indices <- unique(c(spatial_theta_indices, matches))
          }
          if(length(spatial_theta_indices) > 0) {
            # Adjust indices to account for baseline hyperparameters already removed from TETA
            if(is_Surv && (Nsample == 1 || (TRUE %in% sapply(c("rw1", "rw2"), function(x) x %in% unlist(object$basRisk))))) {
              baseline_indices <- grep("baseline", object$misc$theta.tags)
              if(length(baseline_indices) > 0) {
                # Adjust spatial indices by subtracting the number of baseline indices that come before each spatial index
                adjusted_spatial_indices <- sapply(spatial_theta_indices, function(idx) {
                  idx - sum(baseline_indices < idx)
                })
                # Filter out any indices that were baseline (they don't exist in TETA anymore)
                adjusted_spatial_indices <- adjusted_spatial_indices[!spatial_theta_indices %in% baseline_indices]
                spatial_theta_indices <- adjusted_spatial_indices
              }
            }
            # Remove spatial hyperparameters from TETA
            if(length(spatial_theta_indices) > 0) {
              TETA <- TETA[-spatial_theta_indices, , drop=FALSE]
            }
          }
        }
        r <- INLA::inla(formula = formula(FRM3),
                  data = uData,
                  offset = uData$off,
                  E = uData$E..coxph,
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
        r <- .inla_run_many_safe(NsampleHY, wd, num.threads = object$.args$num.threads, cleanup = !TRUE, verbose = !TRUE)
        INLA::inla.setOption(INLAjoint.features=FALSE)
        INLA::inla.setOption(malloc.lib='mi')
        unlink(wd, recursive = TRUE)
        # Only sample IID random effects if there are any left after filtering spatial effects
        if(length(SEL) > 0) {
          RE_values <- do.call(cbind, sapply(1:NsampleHY, function(S) sapply(INLA::inla.posterior.sample(NsampleRE, r[[S]], selection=SEL), function(x) x$latent), simplify=F))
        }else{
          # No IID random effects to sample - all effects are spatial
          RE_values <- NULL
        }
        # Process IID random effects only if they exist
        if(!is.null(RE_values)) {
          if(!object$corLong){
            NRE_i <- length(SEL) # number of random effects
          }else{
            NRE_i <- length(c(object[["REstruc"]], object[["REstrucS"]])) # number of random effects
          }
          NRE_ii <- (dim(RE_values)[1]/NRE_i)/NsampleFE # number of individuals
          id_REV <- data.frame(sapply(1:NsampleFE, function(x) rep(1:NRE_ii, NRE_i) + rep((0:(NRE_i-1))*(NRE_ii*NsampleFE), each=NRE_ii) + rep(NRE_ii*(x-1), (NRE_ii*NRE_i))))
          RE_values <- do.call(cbind, apply(RE_values, 2, function(x) apply(id_REV, 2, function(xx) x[xx]), simplify=F))
          if((NRE_ii + NRE_i)==2 & is_Surv & !is_Long){ # just one vector (may need to adapt for random intercept longitudinal?)
            RE_values <- c(RE_values)
          }
        }
        # Only process IID random effects if they exist
        if(!is.null(RE_values)) {
          # Make grep results robust to handle empty matches
          grep_results <- sapply(c(names_reS, names_reL), function(x) {
            result <- grep(paste0("\\b",x, "\\b"), ct$tag)
            if(length(result) == 0) return(NA) else return(result[1])
          }, USE.NAMES = FALSE)
          if(NRE_ii>1) RE_values <- RE_values[c(sapply(1:(length(unique(ND[,id]))/NsampleFE), function(x) rep(1, NRE_i)+(length(unique(ND[,id]))/NsampleFE)*(seq(1, NRE_i)-1)+(which(unique(ND[,id]) == x)-1))),]
        }
        if(idPredt!=1) idLoopSet <- FALSE else idLoopSet <- TRUE
        if(idLoopSet){ # save all random effects before selecting for each individuals
          if(!is.null(RE_values)) {
            RE_valuesG <- RE_values
          }else{
            # No IID effects sampled via INLA (either no RE or all spatial)
            # Initialize RE_valuesG as an empty matrix to prevent "object not found" errors
            RE_valuesG <- matrix(numeric(0), nrow=0, ncol=0)
          }
        }
        if(!is.null(RE_values) && !is.null(RE_valuesG)){ # Only process if IID effects exist
          Nreps <- trunc(Nsample/Nsample)
          Nadds <- (Nsample %% Nsample)/Nsample
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
      }else if(!exists("RE_valuesG")){
        RE_valuesG <- NULL
      }
      if(is_Long){
        # Handle spatial random effects with set.samples
        RE_valuesSpatial <- NULL
        if(!is.null(set.samples)){
          for(rsmp in 1:length(set.samples)){
            Nrsmp <- names(set.samples)[rsmp]
            if(is.null(dim(set.samples[[rsmp]]))){ # needs to be polished to fit more models
              RE_valuesSpatial <- rep(set.samples[[rsmp]], NsampleRE) # Store spatial values separately
            }else{
              RE_valuesSpatial <- rep(set.samples[[rsmp]][idPredt,], NsampleRE)
            }
          }
          # Only overwrite RE_valuesG if there are NO IID random effects (only spatial)
          if(length(filtered_re_struc) == 0) {
            # Only spatial effects exist - use the original logic
            RE_valuesG <- RE_valuesSpatial
          }
          # Otherwise, keep RE_valuesG as IID values and handle spatial separately
        }
        if(!idLoop){
          # Only process RE_valuesG if it exists and has valid dimensions
          if(!is.null(RE_valuesG) && !is.null(dim(RE_valuesG)) && nrow(RE_valuesG) > 0 && length(unique(ND[, object$id]))/NsampleFE>1){ # only if there are more than 1 individual
            if(!is.null(object[["REstrucS"]])){
              RE_valuesL <- RE_valuesG[-FRAIL_ind,][1:nRE_for_SEL+ rep((RECOUNT_-1)*nRE_for_SEL, nRE_for_SEL),]
            }else{
              RE_valuesL <- RE_valuesG[1:nRE_for_SEL+ rep((RECOUNT_-1)*nRE_for_SEL, nRE_for_SEL),]
            }
          }else{
            RE_valuesL <- RE_valuesG
          }
          ND <- newData[newData[, object$id] == idPred,] # back to individuals now that random effects are done
        }
        if(FEonly && !is.null(RE_values)) RE_values <- matrix(0, nrow = nrow(RE_values), ncol=ncol(RE_values))
        # compute linear predictors for each sample at NtimePoints
        NEWdata[paste0("ID", object[["REstruc"]])] <- NEWdata[paste0("W", object[["REstruc"]])]
        # A matrix for offset computation
        A_LP <- new("dgTMatrix", Dim=c(length(NEWdata[[1]])-length(survPart), sum(ct$length)))
        if(K==1){
          Lout <- 1
        }else{
          Lout <- unique(c(sapply(1:length(object$longOutcome), function(x) grep(paste0(object$longOutcome[[x]], "_L", x), names(NEWdata$Yjoint)))))
        }
        if(is_Surv){
          non_Lout <- setdiff(seq_along(NEWdata$Yjoint), Lout)
          non_Lout_mat <- sapply(NEWdata$Yjoint[non_Lout], function(x) {
            if(inherits(x, "inla.surv")){
              x$time
            }else{
              x
            }
          })
          indL <- which(rowSums(!is.na(non_Lout_mat)) == 0)
        }else if(K>1){
          indL <- c(sapply(NEWdata$Yjoint[Lout], function(x) which(!is.na(x))))
        }else{
          indL <- which(!is.na(NEWdata$Yjoint))
        }
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
          POSc <- rep(1:Nsample, NsampleRE)+rep(0:(NsampleRE-1)*Nsample, each=Nsample)
          # Create RE_mat using IID random effect values
          RE_mat <- new("dgTMatrix",
                        i = as.integer(rep(posRE, length(POSc))-1),
                        j = as.integer(rep(POSc, each=length(posRE))-1),
                        x=c(RE_valuesL), Dim=dim(ParVal2))
          ParVal2 <- ParVal2 + RE_mat
          # Add parallel path for spatial random effects if they exist alongside IID effects
          if(!is.null(set.samples) && !is.null(RE_valuesSpatial) && length(filtered_re_struc) > 0){
            # Get spatial random effect names and positions
            SPATIAL_RE_NAMES <- names(set.samples)
            posSpatialRE <- ct$start[unlist(sapply(SPATIAL_RE_NAMES, function(x) which(ct$tag==x)))]
            if(length(posSpatialRE) > 0) {
              # Add scaling vector of 1s to A_LP for spatial random effects
              A_LP[, posSpatialRE] <- 1
            }
          }
          # ParVal2[posRE, POSc] <- RE_valuesL
          LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal2))
        }else{
          ParVal[posRE, ] <- RE_valuesL
          LP_long <- t(as.matrix(INLA::inla.as.dgTMatrix(A_LP, na.rm=TRUE) %*% ParVal))
        }
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
            famErr_ind <- suppressWarnings(as.integer(gsub(".*\\[(.*)\\].*", "\\1", rownames(object$summary.hyperpar)[hyperr[fer]])))
            if(is.na(famErr_ind)) famErr_ind <- 1
            LP_long[, (1:NTP)+NTP*(famErr_ind-1)] <- LP_long[, (1:NTP)+NTP*(famErr_ind-1)] + REsamF
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
        newPredL <- data.frame(rep(as.factor(rep(idPred, length(LdataPred[, object$id]))), K), rep(LdataPred[, object$timeVar], K),
                               rep(object$longOutcome, each=NTP), RESpredL)
        colnames(newPredL) <- c(object$id, object$timeVar, "Outcome", addNamesL)
        predL <- rbind(predL, newPredL)
        # rm("RE_valuesL")
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
      if(length(survPart)<NTP & length(T_nam)>0){ # only for left truncation? not all time points in survival
        NTP <- length(survPart)
        TPO <- tail(TPO, NTP)
      }
      if(startP==max(TPO) & max(TPO)==horizon) stop(paste0("You ask for predictions at horizon ", horizon, " conditional on data up to time ", startP, ". Please revise horizon or use Csurv argument."))
      TPO2 <- TPO[TPO>=startP]
      NTP2 <- length(TPO2)
      NTP_s <- NTP-NTP2+1
      survPart2 <- survPart[c(unlist(sapply(1:M, function(x) which(NEWdata[[paste0("baseline", x, ".hazard.time")]][NEWdata[[paste0("baseline", x, ".hazard.idx")]]!=0] %in% TPO2))))] # extract part where there is an actual risk
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
            mesh1d <- fmesher::fm_mesh_1d(loc = object$summary.random[[paste0("baseline", m, ".hazard")]]$ID, degree = 1)
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
        # if(length(assocNs2)>0) LP_longs <- LP_long[, -indL][, rep(NTP_s:NTP, length(assocNs2))+rep(patternAsso2-1, each=NTP2)*NTP]
        if(length(assocNa)>0) LP_longs <- LP_long[, -indL][, rep(NTP_s:NTP, length(assocNs2))+rep(patternAsso2-1, each=NTP2)*NTP]
        # I assume all associations are contiguous here (I think it's always true!)
        if(length(assocPos) > 0) {
          ct2$start[assocPos] <- ct2$start[assocPos] - c(0, cumsum(ct2$length[assocPos])[-length(assocPos)]) + cumsum(c(0, rep(NTP2, length(assocNs)-1)))
          ct2$start[-c(1:assocPos[length(assocPos)])] <- ct2$start[-c(1:assocPos[length(assocPos)])] - sum(ct2$length[assocPos]) + NTP2*length(assocNs)#dim(LP_longs)[2]
          ct2$length[assocPos] <- NTP2
        }
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
      if(length(BLpos)>0) A_SP[, c(sapply(BLpos, function(x) ct2$start[x]:(ct2$start[x]+ct2$length[x]-1)))] <- Aproj
      if(exists("assocNa")){
        if(!is.null(assocNa)){          # SET ASSOCIATION INDICATOR HERE INSTEAD OF IS_LONG
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
      }
      if(!is.null(set.samples)){
        for(rsmp in 1:length(set.samples)){
          Nrsmp <- names(set.samples)[rsmp]
          # Check if the spatial effect exists in ct2 before indexing
          if(any(ct2$tag==Nrsmp)) {
            A_SP[, ct2$start[ct2$tag==Nrsmp]] <- 1
          }
          if(is.null(dim(set.samples[[rsmp]]))){
            ParVal[ct$start[ct$tag==Nrsmp], ] <- set.samples[[rsmp]]
          }else{
            ParVal[ct$start[ct$tag==Nrsmp], ] <- set.samples[[rsmp]][idPredt,]
          }
        }
        # need to remove the corresponding random effect from formula
        # if no RE left => skip inla call
        # set sampled values of extra stuff and remove corresponding random effect. (i.e., do not compute posterior)
      }
      # merge
      # scale the association parameters
      if(NsampleRE==F) nsamplere <- 1 else nsamplere <- NsampleRE
      if(NsampleRE!=F & loopRE){
        LP_surv <- NULL
        for(NSRE in 1:NsampleRE){
          # need to fix this since the non loopRE version has been modified!
          # SET ASSOCIATION INDICATOR HERE INSTEAD OF IS_LONG
          if(is_Long && length(assocPos) > 0) {
            SASCP <- t(LP_longs[(1:Nsample + (NSRE-1)*Nsample), ] * sapply(assocNs, function(x) SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==x)])[, rep(1:length(assocPos), each=NTP2)])
            ParValS <- rbind(ParVal[1:(ct$start[assocPos][1]-1), ], SASCP, ParVal[-c(1:(ct$start[assocPos][1] + sum(ct$length[assocPos]) -1)), ])
          }else{
            ParValS <- ParVal
          }
          LP_surv <- rbind(LP_surv, exp(t(as.matrix(INLA::inla.as.dgTMatrix(A_SP, na.rm=TRUE) %*% ParValS))))
        }
      }else{
        # set up matrix to have shared part from longitudinal (LP_longs) scaled by association parameters (assocScaler)
        if(is_Long){
          SASCP <- NULL
          DECAL <- 0 # need to shift when SRE_ind between two time dependent variables
          SRE_indic <- 1
          for(ias in 1:length(assocNs)){
            # Initialize SASCP_t to prevent "object not found" errors
            SASCP_t <- NULL
            if(ias %in% SRE_inda){ # SRE_ind
              assocScaler <- SMPH[, which(gsub("Beta for ", "", colnames(SMPH))==assocNs[ias])][rep(1:Nsample, each=nsamplere)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
              if(!is.null(dim(RE_valuesL))){
                SASCP_t <- RE_valuesL[SRE_inda2[SRE_indic],]*assocScaler # time fixed so only 1 line required
                SRE_indic <- SRE_indic+1
              }else{
                SASCP_t <- RE_valuesL*assocScaler # time fixed so only 1 line required
              }
              PZ <- which(ct2$tag==assocNs[ias]) # position of current assoc
              if(length(PZ) > 0) {
                m_ind <- as.integer(strsplit(assocNs[ias], "_S")[[1]][2])
                PRM <- (ct2$start[PZ]+1):(ct2$start[PZ]+(ct2$length[PZ]-1)) # remove other time points
                if(!(FALSE %in% c(PRM==sort(PRM))) & length(PRM!=2)){
                  A_SP <- A_SP[, -PRM]
                  A_SP[which(!is.na(A_SP[, ct2$start[PZ]])), ct2$start[PZ]] <- 1
                  ct2$start[-c(1:PZ)] <- ct2$start[-c(1:PZ)] - length(PRM)
                  ct2$length[PZ] <- 1
                }
              }
              DECAL <- DECAL + NTP2
            }else{ # other associations
              if(length(which(gsub("Beta for ", "", colnames(SMPH))==assocNs[ias]))>0){
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
          if(length(assocPos) > 0) {
            ParValS <- rbind(ParVal[1:(ct$start[assocPos][1]-1), rep(1:Nsample, each=nsamplere)], SASCP, ParVal[-c(1:(ct$start[assocPos][1] + sum(ct$length[assocPos]) -1)), rep(1:Nsample, each=nsamplere)])
          }else{
            ParValS <- ParVal
          }
        }else{
          ParValS <- ParVal
          PS_nosetup <- TRUE
        }
        # if(reloadCT){ # rerun the code to clean empty columns
        #   ct2_save <- ct2
        #   A_SP_save <- A_SP
        #   reloadCT <- FALSE
        # }
        # ct2 <- ct2_save
        # A_SP <- A_SP_save
        if(!is.null(object[["REstrucS"]])){ # frailty terms?
          # select random effects values for the current individual (only if there are more than 1 individual)
          if(!is.null(RE_valuesG) && !is.null(dim(RE_valuesG)) && nrow(RE_valuesG) > 0 && length(unique(NDS[, object$id]))/NsampleFE>1){ # only if there are more than 1 individual
            RE_valuesS <- RE_valuesG[1:NRE_i+ rep((RECOUNT_-1)*NRE_i, NRE_i),]#[FRAIL_ind,]
            if(!is.null(dim(RE_valuesS))) RE_valuesS <- RE_valuesS[FRAIL_ind,]
          }else if(!is_Long & (length(unique(NDS[, object$id]))/NsampleFE)==1){
            RE_valuesS <- RE_valuesG
          }else{
            # Handle case where RE_valuesG doesn't exist or is invalid
            RE_valuesS <- NULL
          }
          if(exists("PS_nosetup")) ParValS <- ParVal[, rep(1:ncol(ParVal), NsampleRE)]
          if(!is.null(RE_valuesS)) {
            ParValS[ct$start[sapply(object[["REstrucS"]], function(x) which(ct2$tag==x))],] <- RE_valuesS
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
              m_intiCT <- which(ct2$tag == gsub("Beta for ", "", colnames(SMPH)[m_inti]))
              m_ind <- na.omit(sapply(sapply(paste0(object[["REstrucS"]], "_S"), function(x) strsplit(colnames(SMPH)[m_inti], x)[[1]][2]), function(x) as.integer(x)))
              PRM <- (ct2$start[m_intiCT]+1):(ct2$start[m_intiCT]+(ct2$length[m_intiCT]-1)) # remove other time points
              A_SP <- A_SP[, -PRM]
              ParValS <- ParValS[-PRM,]
              A_SP[which(!is.na(A_SP[, ct2$start[m_intiCT]])), ct2$start[m_intiCT]] <- 1
              ct2$start[-c(1:m_intiCT)] <- ct2$start[-c(1:m_intiCT)] - length(PRM)
              ct2$length[m_intiCT] <- 1
              # compute scaled frailty term and insert in shared part
              assocScaler <- SMPH[, m_inti][rep(1:Nsample, nsamplere)]#[rep(1:NTP, M),]*kronecker(assocPoints, matrix(1, ncol=NTP, nrow=NTP))
              # REval_ind <- which(c(object[["REstrucS"]], object[["REstruc"]]) == strsplit(gsub("Beta for ", "", colnames(SMPH)[m_inti]), paste0("_S", m_ind))[[1]])
              if(!is.null(dim(RE_valuesS))){
                ParValS[ct2$start[m_intiCT], ] <- RE_valuesS[ias, ]*assocScaler
              }else{
                ParValS[ct2$start[m_intiCT], ] <- RE_valuesS*assocScaler
              }
            }
          }
        }
        if(exists("Nrsmp")){ # if some samples are set, they may need scaling for shared frailty
          m_inti0 <- which(ct2$tag %in% Nrsmp)
          m_inti1 <- which(ct2$length[m_inti0]>1)
          for(m_intin in m_inti1){ # remove other time points
            PRM <- (ct2$start[m_intin]+1):(ct2$start[m_intin]+(ct2$length[m_intin]-1))
            A_SP <- A_SP[, -PRM]
            ParValS <- ParValS[-PRM,]
            ct2$start[-c(1:m_intin)] <- ct2$start[-c(1:m_intin)] - length(PRM)
            ct2$length[m_intin] <- 1
          }
          if(length(unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))>0){
            for(ias in 1:length(unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH))))))){
              m_inti <- unlist(sapply(gsub("^ID", "", Nrsmp), function(x) grep(paste0(x, "_S"), gsub("Beta for ", "", colnames(SMPH)))))[ias]
              m_intiCT <- which(ct2$tag == gsub("Beta for ", "", colnames(SMPH)[m_inti]))
              m_ind <- na.omit(sapply(sapply(paste0(gsub("^ID", "", Nrsmp), "_S"), function(x) strsplit(colnames(SMPH)[m_inti], x)[[1]][2]), function(x) as.integer(x)))
              PRM <- (ct2$start[m_intiCT]+1):(ct2$start[m_intiCT]+(ct2$length[m_intiCT]-1)) # remove other time points
              A_SP <- A_SP[, -PRM]
              ParValS <- ParValS[-PRM,]
              A_SP[which(!is.na(A_SP[, ct2$start[m_intiCT]])), ct2$start[m_intiCT]] <- 1
              ct2$start[-c(1:m_intiCT)] <- ct2$start[-c(1:m_intiCT)] - length(PRM)
              ct2$length[m_intiCT] <- 1
              # compute scaled frailty term and insert in shared part
              if(is.null(dim(set.samples[[ias]]))){
                ParValS[ct2$start[m_intiCT], ] <- set.samples[[ias]] * SMPH[, m_inti]
              }else{
                ParValS[ct2$start[m_intiCT], ] <- set.samples[[ias]][idPredt,] * SMPH[, m_inti]
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
      newPredS <- data.frame(as.factor(rep(idPred, M*NTP)), rep(TPO, M),
                             rep(paste0("S_", 1:M), each=NTP), Risk2)
      colnames(newPredS) <- c(idname, TimeVar, "Outcome", addNamesS)
      if(survival){
        # compute survival curve
        # take middle of intervals
        SurvSamp <- NULL
        for(m in 1:M){
          rmTP <- c(rep(1:length(TPO2), M-1)+(rep((1:M)[-m], each=length(TPO2))-1)*length(TPO2))
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
        if(length(which(is.nan(SurvSamp[1,])))>0){
          SurvSamp <- SurvSamp[, -which(is.nan(SurvSamp[1,]))]
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
  }
  if(!silentMode) message(paste0("...done!"))
  return(list("PredL"=predL, "PredS"=predS))
}








