#' Fit a multivariate joint model for longitudinal and/or survival data
#'
#' @description This function fits a multivariate joint model for longitudinal and/or survival data.
#' The longitudinal outcomes are modeled with mixed effects models and can handle various distributions.
#' The survival outcomes (i.e., terminal event with possibly competing risks) are modeled with Cox proportional
#' hazards regression models. Various association structures can be specified between the longitudinal and
#' survival outcomes. The inference is based on Integrated Nested Laplace Approximations (Rue et al., 2009).
#'
#' @param formSurv the formula for the time-to-event outcome, with the response given as an inla.surv() object.
#' Keep it as NULL if no survival part is needed and give a list of formulas for competing risks.
#' @param formLong the formula for the longitudinal outcome, structured as with the lme4 package
#' (Mächler et al., 2015) for linear mixed-effects models (i.e., including random effects within parenthesis).
#' Keep it as NULL if no longitudinal part is needed and give a list of formulas for multivariate longitudinal.
#' At the moment, the random effects are limited to intercept, linear slope or random intercept + linear slope.
#' @param dataSurv the dataset for the survival part. Keep it as NULL if no survival part is needed or if
#' the survival data is in the longitudinal dataset (it will extract the last line for each individual as
#' the survival dataset).
#' @param dataLong the dataset for the longitudinal part. Keep it as NULL if no longitudinal part is needed.
#' For multivariate longitudinal models, either give one dataset with all outcomes and covariates or
#' a list of datasets for each longitudinal marker.
#' @param id the name of the variable to identify individuals or grouped repeated measurements.
#' Keep is as NULL if no longitudinal part is needed.
#' @param timeVar a character string (or a vector) giving the name of the time-varying variable(s). Functions of time can
#' be included in formulas, they first need to be set up as a univariate function with name fX,
#' where X is a number between 1 and 20. Then the function can be used directly in the formula (see example below).
#' @param family a character string (or a vector) giving the name of families for the longitudinal outcomes. The
#' list of the available families is given by names(inla.models()$likelihood).
#' @param link a character string (or a vector) giving the link function associated to the families for the longitudinal
#' outcomes. The links available for a family is given in the associated doc: inla.doc("familyName").
#' The link should be a vector of the same size as the family parameter and should be set to "default" for default
#' (i.e., identity for gaussian, log for poisson, logit for bimomial,...).
#' @param basRisk the baseline risk of event (should be a vector in case of competing risks). It can be defined as
#' parametric with either "exponentialsurv" for exponential baseline or "weibullsurv" for Weibull baseline.
#' Alternatively, there are two options to avoid parametric assumptions on the shape of the baseline risk: "rw1"
#' for random walks of order one prior that corresponds to a smooth spline function based
#' on first order differences. The second option "rw2" assigns a random walk order two prior that
#' corresponds to a smooth spline function based on second order differences. This second option
#' provides a smoother spline compared to order one since the smoothing is then done on the second order.
#' @param NbasRisk the number of intervals for the baseline risk function, only one value should be provided
#' and the same number of intervals is used for each risk submodel in cae of competing risks.
#' @param assoc a character string that specifies the association between the longitudinal and survival
#' components. The available options are "CV" for sharing the current value of the linear predictor, "CS"
#' for the current slope, "CV_CS" for the current value and the current slope, "SRE" for shared random effects
#' (i.e., sharing the individual deviation from the mean at time t as defined by the random effects),
#' "SRE_ind" for shared random effect independent (each random effect's individual deviation is associated
#' to an association parameter in the survival submodel) and ""
#' (empty string) for no association. When there are either
#' multiple longitudinal markers or multiple competing events, this should be a vector. In
#' case of both multiple markers and events, it should be a list with one element per longitudinal marker
#' and each element is a vector of association for each competing event. Keep it as NULL to have no
#' association between longitudinal and survival components or if there is no survival component.
#' @param corLong a boolean that only applies when multiple longitudinal markers are fitted: should
#' the random effects  accross markers be correlated (TRUE) or independent (FALSE)? Default is FALSE.
#'
#' @param control a list of control values with components: \describe{
#'
#'   \item{\code{priorFixed}}{list with mean and standard deviations for the Gaussian prior distribution
#'   for the fixed effects. Default is \code{list(mean=0, prec=0.01, mean.intercept=0, prec.intercept=0.01)},
#'   where mean and prec are the mean and precision (i.e., inverse of the variance) of the fixed effects,
#'   respectively and mean.intercept and prec.intercept are the corresponding parameters for the fixed
#'   intercept.}
#'   \item{\code{priorAssoc}}{list with mean and standard deviations for the Gaussian prior distribution
#'   for the association parameters. Default is \code{list(mean=0, prec=0.01)}}
#'   \item{\code{priorRandom}}{list with prior distribution for the multivariate random effects
#'   (Inverse Wishart). Default is \code{list(r=10, R=1)}, see "inla.doc("iidkd") for more details.}
#'   \item{\code{int.strategy}}{a character string giving the strategy for the numerical integration
#'   used to approximate the marginal posterior distributions of the latent field. Available options are
#'   "ccd" (default), "grid" or "eb" (empirical Bayes). The empirical Bayes uses only the mode of the
#'   approximations for the integration, which speed up and simplifies computations. It can be pictured as
#'   a tradeoff between Bayesian and frequentist estimation strategies while the default full Bayesian
#'   accounts for uncertainty by using the mode and the curvature at the mode.}
#'   \item{\code{cpo}}{TRUE/FALSE: Default is FALSE, set to TRUE to compute the Conditional Predictive Ordinate.}
#'   \item{\code{cfg}}{TRUE/FALSE: Default is FALSE, set to TRUE to be able to sample from the posterior
#'   distribution.}
#'   \item{\code{safemode}}{TRUE/FALSE: use the INLA safe mode (automatically reruns in case of negative
#'   eigenvalue(s) in the Hessian, reruns with adjusted starting values in case of crash). Default is TRUE
#'   (activated).
#'   The message ```*** inla.core.safe``` appears when the safe mode is running, it improves the inference
#'   of the hyperparameters and can be ignored. To remove this safe mode, switch the boolean to FALSE (it can
#'   save some computation time but may return slightly less precise estimates for some hyperparameters).
#'   }
#'   \item{\code{verbose}}{TRUE/FALSE: prints details of the INLA algorithm. Default is FALSE.}
#'   \item{\code{keep}}{TRUE/FALSE: keep internal files. Default is FALSE. (expert option)}
#'}
#'
#'
#'
#' @return An object of class \code{INLAjoint}. See \code{\link{INLAjoint.object}} for
#'   details.
#' @references
#' Rustand, D., van Niekerk, J., Teixeira Krainski, E., Rue, H. and Proust-Lima, C. (2022).
#' Fast and flexible inference approach for joint models of multivariate longitudinal and
#' survival data using Integrated Nested Laplace Approximations.
#' https://arxiv.org/abs/2203.06256
#'
#' Rustand, D., van Niekerk, J., Rue, H., Tournigand, C., Rondeau, V. and Briollais, L. (2021).
#' Bayesian Estimation of Two-Part Joint Models for a Longitudinal Semicontinuous Biomarker
#' and a Terminal Event with R-INLA: Interests for Cancer Clinical Trial Evaluation.
#' https://arxiv.org/abs/2010.13704
#'
#' Rue, H., Martino, S. and Chopin, N. (2009). Approximate Bayesian inference for latent
#' Gaussian models by using integrated nested Laplace approximations. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology), 71: 319-392.
#' https://doi.org/10.1111/j.1467-9868.2008.00700.x
#'
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects
#' Models Using lme4. Journal of Statistical Software, 67(1), 1–48.
#' https://doi.org/10.18637/jss.v067.i01
#'
#' Contact: \email{INLAjoint@gmail.com}
#'@export
#'
#' @examples
#'
#' # joint model with 3 longitudinal / 3 competing risks of event
#' data(Long)
#' data(Surv)
#' YD1 <- inla.surv(time = c(Surv$deathTimes), event = c(Surv$Event1)) # Event 1
#' YD2 <- inla.surv(time = c(Surv$deathTimes), event = c(Surv$Event2)) # Event 2
#' YD3 <- inla.surv(time = c(Surv$deathTimes), event = c(Surv$Event3)) # Event 3
#' f1 <- function(x) x^2 # quadratic function of time for first marker
#' Nsplines <- ns(Long$time, knots=2) # 2 ns splines for second marker
#' f2 <- function(x) predict(Nsplines, x)[,1]
#' f3 <- function(x) predict(Nsplines, x)[,2]
#'
#' JMINLA <- joint(
#'  formLong = list(Y1 ~ time + f1(time) + ctsX + binX + (1 + time + f1(time) | Id),
#'                  Y2 ~ time + f2(time) + f3(time) + ctsX + binX + (1 | Id),
#'                  Y3 ~ time + ctsX + binX + (1 | Id)), dataLong = Long,
#'  formSurv = list(YD1 ~ binX + ctsX,
#'                  YD2 ~ binX,
#'                  YD3 ~ ctsX),
#'  id = "Id", timeVar = "time", corLong=TRUE,
#'  family = c("gaussian", "poisson", "binomial"), basRisk = c("rw1", "rw1", "rw1"),
#'  assoc = list(c("CV", "CS", ""),  c("CV", "", "SRE"), c("", "CV", "")),
#'  control=list(int.strategy="eb"))
#'
#' summary(JMINLA)
#' # 'sdcor' to switch from variance-covariance to standard
#' # deviation-correlation and 'hazr' to switch parameters
#' # in survival submodels from mean to hazard ratios (exp(mean)).
#' summary(JMINLA, sdcor=TRUE, hazr=TRUE)
#'

joint <- function(formSurv = NULL, formLong = NULL, dataSurv=NULL, dataLong=NULL,
                  id=NULL, timeVar=NULL, family = "gaussian", link = "default",
                  basRisk = "rw1", NbasRisk = 15, assoc = NULL, corLong=FALSE,
                  control = list(), ...) {

  is_Long <- !is.null(formLong) # longitudinal component?
  is_Surv <- !is.null(formSurv) # survival component?

  # Number of survival events = M and conversion to list if M=1
  if(is_Surv){
    if (!is.list(formSurv)) {
      if(class(formSurv)!="formula") stop("formSurv must be a formula or a list of formulas")
      formSurv <- list(formSurv)
      M <- 1
    }else{
      M <- length(formSurv) # number of time-to-event outcomes
      for(m in 1:M){
        if(class(formSurv[[m]])!="formula") stop("formSurv must be a formula or a list of formulas")
      }
    }
    if(length(basRisk)!=M) stop(paste0("basrisk must contain a vector of elements with the baseline risk function for
                                       each survival component (i.e., ",M," components while I found ",length(basRisk),
                                       " component(s)."))
    # Check length of baseline risk and conversion to a list
    if(length(basRisk)==1 & !is.list(basRisk)){
      basRisk <- list(basRisk)
    }else if(length(basRisk)>1 & !is.list(basRisk)){
      basRisk <- as.list(basRisk)
    }else if(length(basRisk)>1 & is.list(basRisk)){
      if(length(basRisk) != M) stop(paste0("The length of basRisk must match the number of formulas for the survival part (found  ", length(basRisk), " items in basRisk and ", M, " formulas for survival)."))
    }
  }
  if(is_Long){
    if(!is.list(formLong)){ # Number of longitudinal markers = K and conversion to list if K=1
      if(class(formLong)!="formula") stop("formLong must be a formula or a list of formulas")
      formLong <- list(formLong)
      K <- 1
    }else{
      K <- length(formLong) # number of markers
      for(k in 1:K){
        if(class(formLong[[k]])!="formula") stop("formLong must be a formula or a list of formulas")
      }
    }
    # Check length of family and conversion to list
    if (!is.list(family)) {
      family <- as.list(family)
      if(length(family)!=K) stop(paste("number of families must match number of longitudinal markers (i.e., ", K, ")"))
      if(!length(link)==1){
        if(length(family)!=length(link)) stop(paste("number of families must match number of link functions (found ", length(family), " families and ", length(link), "link functions."))
      }else if(K>1){
        link <- rep("default", K) # if link not provided, assume default
      }
    }
    # Either one dataset per formula or one dataset for all markers
    if(class(dataLong)=="list"){
      if(length(dataLong)%in%c(K, 1)){
        if(length(dataLong)==1) oneData=TRUE else oneData=FALSE # only one Dataset for longitudinal
      }else{
        stop(paste("The number of dataset for the longitudinal markers must be either one or equal to the number of markers (i.e., ", K, ")"))
      }
    }else{
      oneData=TRUE
      dataLong <- list(dataLong)
    }

    # check timeVar
    if(length(timeVar)>1) stop("timeVar must only contain the time variable name.")

    if(is_Long){
      # replace special characters in factor variables
      for(i in 1:length(dataLong)){
        colClass <- sapply(dataLong[[i]], class)
        #dataLong[[i]][,which(colClass=="character")] <- sapply(dataLong[[i]][,which(colClass=="character")], function(x) sub("-","", x))
        dataLong[[i]][,which(colClass=="character")] <- sapply(dataLong[[i]][,which(colClass=="character")], function(x) sub("[^[:alnum:] ]","", x))
        if(length(which(colClass=="factor"))>0){
          for(fctrs in 1:length(which(colClass=="factor"))){
            lvlFact <- levels(dataLong[[i]][,which(colClass=="factor")[fctrs]]) # save reference level because otherwise it can change it
            #dataLong[[i]][,which(colClass=="factor")[fctrs]] <- factor(sub("-","", dataLong[[i]][,which(colClass=="factor")[fctrs]]), levels=sub("-","", lvlFact))
            dataLong[[i]][,which(colClass=="factor")[fctrs]] <- factor(sub("[^[:alnum:] ]","", dataLong[[i]][,which(colClass=="factor")[fctrs]]), levels=sub("[^[:alnum:] ]","", lvlFact))
          }
        }
      }
    }

    dataL <- dataLong[[1]] # dataL contains the dataset for marker k (always the same if only one dataset provided)
    if(is.null(timeVar)) warning("There is no time variable in the longitudinal model? (timeVar argument)")
    if(is.null(id)){
      warning("There is no id variable in the longitudinal model? (id argument)")
    }else{
      if(max(as.integer(dataL[,id]))!=length(unique(dataL[,id]))){ # avoid missing ids
        dataL[,id] <- as.integer(as.factor(dataL[,id]))
      }
    }
  }else if(!is_Surv){
    stop("Error: no longitudinal or survival part detected...")
  }

  if(is_Surv){
    # get a dataset with unique line for each ID in case some covariates from the longitudinal
    # part for the association are not provided in the survival model
    if(length(dataSurv)==0){ # if dataSurv not provided, extract it from dataLong
      oneDataS <- TRUE
      LSurvdat <- dataL[c(which(diff(as.numeric(dataL[,which(colnames(dataL)==id)]))==1),
                          length(dataL[,which(colnames(dataL)==id)])),]
      dataSurv <- list(LSurvdat)
    }else if(length(dataSurv)>0){
      # make data as a list
      if(!class(dataSurv)=="list") dataSurv <- list(dataSurv)
      # indicator for one unique survival dataset vs one dataset per model
      if(length(dataSurv)==1) oneDataS <- TRUE else oneDataS <- FALSE
      if(exists("dataL")) LSurvdat <- dataL[c(which(diff(as.numeric(dataL[,which(colnames(dataL)==id)]))==1),
                                              length(dataL[,which(colnames(dataL)==id)])),]
    }
    for(i in 1:length(dataSurv)){
      colClass <- sapply(dataSurv[[i]], class)
      #dataSurv[[i]][,which(colClass=="character")] <- sapply(dataSurv[[i]][,which(colClass=="character")], function(x) sub("-","", x))
      dataSurv[[i]][,which(colClass=="character")] <- sapply(dataSurv[[i]][,which(colClass=="character")], function(x) sub("[^[:alnum:] ]","", x))
      if(length(which(colClass=="factor"))>0){
        for(fctrs in 1:length(which(colClass=="factor"))){
          lvlFact <- levels(dataSurv[[i]][,which(colClass=="factor")[fctrs]]) # save reference level because otherwise it can change it
          #dataSurv[[i]][,which(colClass=="factor")[fctrs]] <- factor(sub("-","", dataSurv[[i]][,which(colClass=="factor")[fctrs]]), levels=sub("-","", lvlFact))
          dataSurv[[i]][,which(colClass=="factor")[fctrs]] <- factor(sub("[^[:alnum:] ]","", dataSurv[[i]][,which(colClass=="factor")[fctrs]]), levels=sub("[^[:alnum:] ]","", lvlFact))        }
      }
      if(!is.null(id)){
        if(id %in% colnames(dataSurv[[i]])){
          if(max(as.integer(dataSurv[[i]][,id]))!=length(unique(dataSurv[[i]][,id]))){ # avoid missing ids
            dataSurv[[i]][,id] <- as.integer(as.factor(dataSurv[[i]][,id]))
          }
        }
      }
    }
    if(!exists("LSurvdat")) LSurvdat <- dataSurv[[1]]
  }


  # Check if no survival => no assoc
  if(!is_Surv & length(assoc)!=0) stop("There is no survival component, therefore assoc should be set to NULL.")
  # Check association structure length and values and conversion to a list
  # It should be a list of K vectors of size M
  # or a vector of size K if M=1
  # or a vector of size M if K=1
  if(length(assoc)!=0){
    if (!is.list(assoc)){
      if(M>1){
        assoc <- list(assoc)
      }else if(M==1){
        assoc <- as.list(assoc)
      }
    }
    if(length(assoc)!=K) stop("Please provide one association per longitudinal marker")
    for(k in 1:K){
      if(length(assoc[[k]])!=M) stop(paste0("Please provide one association per survival outcome for marker ", k))
      for(m in 1:M){
        if(!(assoc[[k]][m] %in% c("CV","CS", "CV_CS", "SRE", "SRE_ind", ""))) stop(paste0("Please choose among the available association structures (NULL, CV, CS, CV_CS, SRE, SRE_ind) for marker ", k, " and event ", m))
      }
    }
  }

  if(!(corLong %in% c(T, F))) stop("corLong must be either TRUE of FALSE")
  # control variables
  if(is.null(control$priorFixed$mean)) control$priorFixed$mean <- 0
  if(is.null(control$priorFixed$prec)) control$priorFixed$prec <- 0.01
  if(is.null(control$priorFixed$mean.intercept)) control$priorFixed$mean.intercept <- 0
  if(is.null(control$priorFixed$prec.intercept))  control$priorFixed$prec.intercept <- 0.01
  if(is.null(control$priorAssoc$mean)) control$priorAssoc$mean <- 0
  if(is.null(control$priorAssoc$prec)) control$priorAssoc$prec <- 0.01
  if(is.null(control$PriorRandom$r)) control$PriorRandom$r <- 10
  if(is.null(control$PriorRandom$R)) control$PriorRandom$R <- 1

  safemode <- ifelse("safemode" %in% names(control), control$safemode, T)
  verbose <- ifelse("verbose" %in% names(control), control$verbose, F)
  cpo <- ifelse("cpo" %in% names(control), control$cpo, F)
  keep <- ifelse("keep" %in% names(control), control$keep, F)
  variant <- ifelse("variant" %in% names(control), control$variant, 1) # for weibull baseline hazard


  int.strategy <- ifelse("int.strategy" %in% names(control), control$int.strategy, "ccd")
  cfg <- ifelse("cfg" %in% names(control), control$cfg, FALSE)
  cpo <- ifelse("cpo" %in% names(control), control$cpo, FALSE)
  assocInit<- 0.1 # switch from default 1 to 0.1 for more stability

  NFT <- 20 # maximum number of functions of time (f1, f2, ...)
  ################################################################# survival part
  cureVar <- NULL
  if(is_Surv){
    modelYS <- vector("list", M) # models for survival outcomes
    data_cox <- vector("list", M) # data for survival outcomes + association terms
    re.weight <- vector("list", M) # data for survival outcomes + association terms
    id_cox <- vector("list", M) # data for survival outcomes + association terms
    ns_cox <- vector("list", M) # data for survival outcomes + association terms
    formAddS <- vector("list", M) # store formula part for random effects in survival if any
    cureVar <- vector("list", M) # store mixture cure predictors names
    IDres <- 0
    REstrucS=NULL # used to have the structure of random effects for survival in output
    if(is_Long) IDassoc <- vector("list", K) # unique identifier for the association between longitudinal and survival
    if(oneDataS) dataS <- dataSurv[[1]]
    IDas <- 0 # to keep track of unique id for association
    for(m in 1:M){ # loop over M survival outcomes
      if(corLong != TRUE) IDres <- 0 # to keep track of unique id for random effects
      if(!oneDataS | m==1){# remove special character "-" from variables modalities
        colClass <- sapply(dataSurv[[m]], class)
        dataSurv[[m]][,which(colClass=="character")] <- sapply(dataSurv[[m]][,which(colClass=="character")], function(x) sub("-","", x))
        if(length(which(colClass=="factor"))>0){
          for(fctrs in 1:length(which(colClass=="factor"))){
            lvlFact <- levels(dataSurv[[m]][,which(colClass=="factor")[fctrs]]) # save reference level because otherwise it can change it
            dataSurv[[m]][,which(colClass=="factor")[fctrs]] <- factor(sub("-","", dataSurv[[m]][,which(colClass=="factor")[fctrs]]), levels=sub("-","", lvlFact))
          }
        }
      }
      if(!oneDataS | m==1) dataS <- dataSurv[[m]]
      # first set up the data and formula for marker m
      modelYS[[m]] <- setup_S_model(formSurv[[m]], formLong, dataS, LSurvdat, timeVar, assoc, id, m, K, M, NFT, corLong)
      # then do the cox expansion to have intervals over the follow-up,
      # these intervals have 2 use: the evaluation of the Bayesian smoothing splines for the baseline risk
      # account for time-dependent component in the association parameters
      # if(length(grep("Intercept", modelYS[[m]][[2]]))!=0){
      #   cstr=TRUE
      # }else{
        cstr=FALSE
      # }
      if(basRisk[[m]]%in%c("rw1", "rw2") & !is.null(modelYS[[m]][[1]][[1]]$cure)){
        # not using warning() function because we want this message to by systematically printed
        print("Warning: Mixture cure model not available for random walks 1 and 2 baseline risk. Please switch to parametric baseline (i.e., exponentialsurv or weibullsurv) to enable mixture cure.")
        modelYS[[m]][[1]][[1]]$cure = NULL
      }
      if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv") & !is.null(modelYS[[m]][[1]][[1]]$cure) & is.null(colnames(modelYS[[m]][[1]][[1]]$cure))){
        # not using warning() function because we want this message to by systematically printed
        print("Warning: Variables names in mixture cure regression model not given, automatically assigning names ('Cure1', 'Cure2', etc.)")
        colnames(modelYS[[m]][[1]][[1]]$cure) <- paste0("Cure", 1:dim(modelYS[[m]][[1]][[1]]$cure)[2], "_S",m)
      }else if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv") & !is.null(modelYS[[m]][[1]][[1]]$cure) & !is.null(colnames(modelYS[[m]][[1]][[1]]$cure))){
        colnames(modelYS[[m]][[1]][[1]]$cure) <- paste0(colnames(modelYS[[m]][[1]][[1]]$cure), "(cure)", "_S",m)
      }
      if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv") & !is.null(modelYS[[m]][[1]][[1]]$cure)) cureVar[[m]] <- colnames(modelYS[[m]][[1]][[1]]$cure)
      if(!is.null(assoc)) YS_assoc <- unlist(assoc[1:K])[seq(m, K*M, by=M)] else YS_assoc <- NULL # extract K association terms associated to time-to-event m
      if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv") & !TRUE %in% c(c("CV", "CS", "CV_CS", "SRE") %in% YS_assoc)){
        DatParam <- data.frame(modelYS[[m]][[1]][[1]]$event, ifelse(modelYS[[m]][[1]][[1]]$time>modelYS[[m]][[1]][[1]]$lower, modelYS[[m]][[1]][[1]]$time, modelYS[[m]][[1]][[1]]$lower),modelYS[[m]][[1]][-1])
        colnames(DatParam)[1] <- paste0("y", m, "..coxph")
        colnames(DatParam)[2] <- paste0("surv", m, "time")
        assign(paste0("cox_event_", m), list(formula=modelYS[[m]][[2]], data=DatParam))
        if(!is.null(modelYS[[m]][[1]][[1]]$cure)) assign(paste0("cure_", m), modelYS[[m]][[1]][[1]]$cure) else assign(paste0("cure_", m), NULL)
        assign(paste0("formS", m), get(paste0("cox_event_", m))$formula)
      }else{
        BR=ifelse(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv"), "rw1", basRisk[[m]])
        assign(paste0("cox_event_", m), # dynamic name to have a different object for each survival outcome m
               inla.coxph(modelYS[[m]][[2]], control.hazard=list(
                 model=BR, scale.model=TRUE,
                 diagonal=1e-2,constr=cstr, n.intervals=NbasRisk,
                 hyper=list(prec=list(prior='pc.prec', param=c(0.5,0.01), initial=3))),
                 data = modelYS[[m]][[1]], tag=as.character(m)))
        if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv")){
          assign(paste0("formS", m), formula(paste0("Yjoint ~", strsplit(as.character(get(paste0("cox_event_",m))$formula)[[3]], paste0("\\+ f\\(baseline", m))[[1]][1], "-1")))
          if(!is.null(modelYS[[m]][[1]][[1]]$cure)) assign(paste0("cure_", m), modelYS[[m]][[1]][[1]]$cure) else assign(paste0("cure_", m), NULL)
        }else{
          assign(paste0("formS", m), get(paste0("cox_event_", m))$formula)
          assign(paste0("cure_", m), NULL)
        }
      }
      ns_cox[[m]] = dim(get(paste0("cox_event_", m))$data)[1] # size of survival part after decomposition of time into intervals
      if(is.null(id_cox[[m]])) id_cox[[m]] <- as.integer(unname(unlist(get(paste0("cox_event_", m))$data[length(get(paste0("cox_event_", m))$data)]))) # repeated individual id after cox expansion
      # weight for time dependent components = middle of the time interval
      if(basRisk[[m]]%in%c("exponentialsurv", "weibullsurv") & !TRUE %in% c(c("CV", "CS", "CV_CS", "SRE") %in% YS_assoc)){ # set event time as weight if parametric, otherwise use middle of interval
        re.weight[[m]] <- ifelse(modelYS[[m]][[1]][[1]]$time>modelYS[[m]][[1]][[1]]$lower, modelYS[[m]][[1]][[1]]$time, modelYS[[m]][[1]][[1]]$lower)
      }else{
        re.weight[[m]] <- unname(unlist(get(paste0("cox_event_", m))$data[paste0("baseline", m, ".hazard.time")] + 0.5 *get(paste0("cox_event_", m))$data[paste0("baseline", m, ".hazard.length")]))
      }
      # set up unique id for association
      if(length(assoc)!=0){
        if(IDas==0 & m==1){ # do this only once
          for(k in 1:K){ # store unique id for each association term (must be unique at each time point instead of individual repeated id)
            # for each marker k, each association term (CV, CS or SRE) must have an unique id for each row of data_cox
            if("CV_CS" %in% assoc[[k]]){
              if(!("CV" %in% assoc[[k]])){
                IDassoc[[k]] <- append(IDassoc[[k]], list("CV"=(IDas + 1:ns_cox[[m]])))
                IDas <- IDas+ns_cox[[m]]
              }
              if(!("CS" %in% assoc[[k]])){
                IDassoc[[k]] <- append(IDassoc[[k]], list("CS"=(IDas + 1:ns_cox[[m]])))
                IDas <- IDas+ns_cox[[m]]
              }
            }
            if("CV" %in% assoc[[k]]){
              IDassoc[[k]] <- append(IDassoc[[k]], list("CV"=(IDas + 1:ns_cox[[m]])))
              IDas <- IDas+ns_cox[[m]]
            }
            if("CS" %in% assoc[[k]]){
              IDassoc[[k]] <- append(IDassoc[[k]], list("CS"=(IDas + 1:ns_cox[[m]])))
              IDas <- IDas+ns_cox[[m]]
            }
            if("SRE" %in% assoc[[k]]){
              IDassoc[[k]] <- append(IDassoc[[k]], list("SRE"=(IDas + 1:ns_cox[[m]])))
              IDas <- IDas+ns_cox[[m]]
            } # SRE_ind is set up with 'copy', which doesn't require an unique ID.
          }
        }
        YS_assoc <- unlist(assoc[1:K])[seq(m, K*M, by=M)] # extract K association terms associated to time-to-event m
        data_cox[[m]] <- get(paste0("cox_event_", m))$data # store the data in this object, it is easier to manipulate compared to object with dynamic name
        for(k in 1:length(YS_assoc)){ # update association id to make them unique instead of individually repeated, so we can account for time dependency
          if(YS_assoc[k]%in%c("CV", "CS", "SRE")){ # one vector
            if(dim(data_cox[[m]][paste0(YS_assoc[k], "_L", k, "_S", m)])[[1]] == length(unlist(IDassoc[[k]][YS_assoc[k]]))){
              data_cox[[m]][paste0(YS_assoc[k], "_L", k, "_S", m)] <- IDassoc[[k]][YS_assoc[k]]
            }else{
              if("CV" == YS_assoc[k]){
                IDassoc[[k]]$CV <- IDas + 1:ns_cox[[m]]
                IDas <- IDas+ns_cox[[m]]
              }
              if("CS" == YS_assoc[k]){
                IDassoc[[k]]$CS <- IDas + 1:ns_cox[[m]]
                IDas <- IDas+ns_cox[[m]]
              }
              if("SRE" == YS_assoc[k]){
                IDassoc[[k]]$SRE <- IDas + 1:ns_cox[[m]]
                IDas <- IDas+ns_cox[[m]]
              }
              data_cox[[m]][paste0(YS_assoc[k], "_L", k, "_S", m)] <- IDassoc[[k]][YS_assoc[k]]
            }
          }else if(YS_assoc[k]=="CV_CS"){ # two vectors for current value and current slope
            if(dim(data_cox[[m]][paste0("CV_L", k, "_S", m)])[[1]] == length(unlist(IDassoc[[k]]["CV"]))){
              data_cox[[m]][paste0("CV_L", k, "_S", m)] <- IDassoc[[k]]["CV"]
            }else{
              IDassoc[[k]]$CV <- IDas + 1:ns_cox[[m]]
              IDas <- IDas+ns_cox[[m]]
              data_cox[[m]][paste0("CV_L", k, "_S", m)] <- IDassoc[[k]]["CV"]
            }
            if(dim(data_cox[[m]][paste0("CS_L", k, "_S", m)])[[1]] == length(unlist(IDassoc[[k]]["CS"]))){
              data_cox[[m]][paste0("CS_L", k, "_S", m)] <- IDassoc[[k]]["CS"]
            }else{
              IDassoc[[k]]$CS <- IDas + 1:ns_cox[[m]]
              IDas <- IDas+ns_cox[[m]]
              data_cox[[m]][paste0("CS_L", k, "_S", m)] <- IDassoc[[k]]["CS"]
            }
          }
        }
      }else{
        data_cox[[m]] <- get(paste0("cox_event_", m))$data # store the data in this object, it is easier to manipulate compared to object with dynamic name
      }
      if(!is.null(modelYS[[m]]$RE_matS)){ # random effects in survival model m
        for(j in 1:ncol(modelYS[[m]]$RE_matS)){
          if(!(colnames(modelYS[[m]]$RE_matS)[j] %in% c(timeVar, c(paste0("f", 1:NFT, timeVar))))){
            if(colnames(modelYS[[m]]$RE_matS)[j]=="Intercept"){
              assign(paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), id_cox[[m]]) # assign variable with dynamic name for random effect
              assign(paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), rep(1, length(id_cox[[m]]))) # assign variable with dynamic name for associated weight
              lid <- length(unique(id_cox[[m]])) # length id
            }else{
              if(exists("data_cox")){
                idVar <- unname(sapply(modelYS[[m]][[2]][,j], function(x) correspondID[which(correspondID[,1]==x),2])) #set id for random effect
                assign(paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), idVar) # assign variable with dynamic name for random effect
                assign(paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), c(rep(1, length(idVar)))) # assign variable with dynamic name for associated weight
                lid <- length(unique(idVar)) # length id
              }else{
                correspondID <- cbind(unique(modelYS[[m]][[2]][,j]), 1:length(unique(modelYS[[m]][[2]][,j])))
                idVar <- unname(sapply(modelYS[[m]][[2]][,j], function(x) correspondID[which(correspondID[,1]==x),2])) #set id for random effect
              }
            }
          }else{
            assign(paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), id_cox[[m]]) # assign variable with dynamic name for random effect
            assign(paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m), c(unname(modelYS[[m]][[2]][,j]))) # assign variable with dynamic name for associated weight
          }
          data_cox[[m]] <- cbind(data_cox[[m]], get(paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m)))
          names(data_cox[[m]]) <- c(names(data_cox[[m]])[-length(names(data_cox[[m]]))], paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m))
          data_cox[[m]] <- cbind(data_cox[[m]], get(paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m)))
          names(data_cox[[m]]) <- c(names(data_cox[[m]])[-length(names(data_cox[[m]]))], paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m))
          if(j==1){
            formAddS[[m]] <- paste0("Yjoint ~ . + ", paste("f(", paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m),",", paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m),", model = 'iid',
                n =", lid,", hyper=list(prec=list(prior='loggamma', param=c(0.01,0.01))))"))
          }else{
            formAddS[[m]] <- update(formAddS[[m]], paste0("Yjoint ~ . + ", paste("f(", paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m),",", paste0("W",colnames(modelYS[[m]]$RE_matS)[j], "_S",m),", model = 'iid',
                n =", lid,", hyper=list(prec=list(prior='loggamma', param=c(0.01,0.01))))")))
          }
        }
        REstrucS <- c(REstrucS, paste0("ID",colnames(modelYS[[m]]$RE_matS)[j], "_S",m))
      }
    }
    dlCox <- NULL # just need to grab the "data.list" from the cox_event object with dynamic name in order to merge it with the rest of the data
    if(M>1){
      dlCoxtemp <- mget(paste0("cox_event_", 1:M))
      for(m in 1:M){
        dlCox <- append(dlCox, dlCoxtemp[[m]][["data.list"]])
      }
    }else{
      dlCox <- append(dlCox, cox_event_1[["data.list"]])
    }
  }
  ################################################################# longitudinal part
  if(is_Long){
    modelYL <- vector("list", K) # outcomes list
    modelFE <- vector("list", K) # fixed effects for each marker k (list of size K)
    modelRE <- vector("list", K) # random effects for each marker k (list of size K)
    Nid <- vector("list", K) # number of individuals for each k
    YL <- list() # to store the outcomes
    dataFE <- NULL # data for the fixed effects part of the longitudinal submodels
    dataRE <- NULL # data for the random effects part of the longitudinal submodels
    dataYL <- NULL # data for the longitudinal outcomes
    Vasso <- NULL
    NAvect <- 0 # vector of NA to fill the vectors up until marker k's part (length of k-1 markers + association)
    IDre <- 0
    fam <- NULL # set up families
    famCtrl <- NULL
    for(k in 1:K){
      if(corLong != TRUE) IDre <- 0 # to keep track of unique id for random effects
      if(!oneData | k==1){# remove special character "-" from factors/character variables modalities
        colClass <- sapply(dataLong[[k]], class)
        dataLong[[k]][,which(colClass=="character")] <- sapply(dataLong[[k]][,which(colClass=="character")], function(x) sub("-","", x))
        if(length(which(colClass=="factor"))>0){
          for(fctrs in 1:length(which(colClass=="factor"))){
            lvlFact <- levels(dataLong[[k]][,which(colClass=="factor")[fctrs]]) # save reference level because otherwise it can change it
            dataLong[[k]][,which(colClass=="factor")[fctrs]] <- factor(sub("-","", dataLong[[k]][,which(colClass=="factor")[fctrs]]), levels=sub("-","", lvlFact))
          }
        }
      }
      if(!oneData | k==1) dataL <- dataLong[[k]] # dataL contains the dataset for marker k (always the same if only one dataset provided)
      if(max(as.integer(dataL[,id]))!=length(unique(dataL[,id]))){ # avoid missing ids
        dataL[,id] <- as.integer(as.factor(dataL[,id]))
      }
      modelYL[[k]] <- setup_Y_model(formLong[[k]], dataL, family[[k]], k) # prepare outcome part for marker k
      modelFE[[k]] <- setup_FE_model(formLong[[k]], dataL, timeVar, k) # prepare fixed effects part for marker k
      modelRE[[k]] <- setup_RE_model(formLong[[k]], dataL, k) # prepare random effects part for marker k
      Nid[[k]] <- length(unique(as.integer(dataL[,id]))) # number of individuals for marker k
      for(j in 1:length(modelFE[[k]][[1]])){ # fixed effects
        if(length(assoc)!=0){ # set up the association for marker k
          Vasso <- NULL # vector of all the association parts
          assoInfo <- NULL
          assoCur <- NULL # current association
          for(m in 1:M){
            assoCur <- assoc[[k]][m]
            if(!dim(data_cox[[m]])[1] %in% assoInfo[which(assoCur == assoInfo[,1]),2]){
              if(assoCur=="CV_CS"){ # if one of CV or CS is already done, don't redo it!
                if(dim(data_cox[[m]])[1] %in% assoInfo[which("CV" == assoInfo[,1]),2]) assoCur <- "CS"
                if(dim(data_cox[[m]])[1] %in% assoInfo[which("CS" == assoInfo[,1]),2]) assoCur <- "CV"
              }
              if(!(TRUE %in% (unlist(strsplit(modelFE[[k]][[1]][j], ".X.")) %in% c(timeVar, c(paste0("f", 1:NFT, timeVar)))))){ # if variable j does not contain a  time-dependent variable
                if(modelFE[[k]][[1]][j]=="Intercept"){ # if intercept
                  if(assoCur %in% c("CV", "CV_CS")){
                    Vasso <- c(Vasso, rep(1, ns_cox[[m]]))
                  }
                }else{
                  # else put corresponding value of the fixed effect variable for the association part
                  if(assoCur %in% c("CV", "CV_CS")){
                    Vasso <- c(Vasso, data_cox[[m]][, modelFE[[k]][[1]][j]])
                  }
                }
                if(assoCur %in% c("CS", "CV_CS", "SRE")){
                  Vasso <- c(Vasso, rep(NA, ns_cox[[m]]))
                }
              }else{ # if time varying variable is in variable j (either alone or with interaction)
                if("X" %in% unlist(strsplit(modelFE[[k]][[1]][j], "\\."))){ # if time variable has an interaction
                  # first identify the time variable
                  if(timeVar %in% (unlist(strsplit(modelFE[[k]][[1]][j], "\\.")))){
                    tvar <- which(unlist(strsplit(modelFE[[k]][[1]][j], "\\.")) %in% timeVar) # useless line?
                    # then identify the other non time-varying variable(s) of the interaction
                    ntvar <-  unlist(strsplit(modelFE[[k]][[1]][j], "\\."))[-which(unlist(strsplit(modelFE[[k]][[1]][j], "\\.")) %in% c(timeVar, "X"))]
                    if(assoCur %in% c("CV", "CV_CS")){
                      Vasso <- c(Vasso, re.weight[[m]] * data_cox[[m]][, ntvar])
                    }
                    if(assoCur %in% c("CS", "CV_CS")){
                      Vasso <- c(Vasso, data_cox[[m]][, ntvar])
                    }
                  }else{
                    # in case of interaction of a function of time, first identify the function of time position in the interaction
                    tvar <- unlist(strsplit(modelFE[[k]][[1]][j], "\\."))[which(unlist(strsplit(modelFE[[k]][[1]][j], "\\.")) %in% c(paste0("f", 1:NFT, timeVar)))]
                    # then identify the other non time-varying variable(s) of the interaction
                    ntvar <-  unlist(strsplit(modelFE[[k]][[1]][j], "\\."))[-which(unlist(strsplit(modelFE[[k]][[1]][j], "\\.")) %in% c(paste0("f", 1:NFT, timeVar),timeVar, "X"))]
                    if(assoCur %in% c("CV", "CV_CS")){
                      # evaluate f function of time at time points re.weight for current value association in survival
                      # and multiply by the other variable for the interaction
                      Vasso <- c(Vasso, unname(sapply(re.weight[[m]], paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == tvar))))*data_cox[[m]][, ntvar])
                    }
                    if(assoCur %in% c("CS", "CV_CS")){
                      # evaluate derivative of f function of time at time points re.weight for current value association in survival
                      DerivValue <- numDeriv::grad(get(paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == tvar))), re.weight[[m]])
                      # and multiply by the other variable for the interaction
                      Vasso <- c(Vasso, DerivValue * data_cox[[m]][, ntvar])
                    }
                  }
                }else{
                  if(modelFE[[k]][[1]][j] == timeVar){
                    if(assoCur %in% c("CV", "CV_CS")){
                      Vasso <- c(Vasso, re.weight[[m]])
                    }
                    if(assoCur %in% c("CS", "CV_CS")){
                      Vasso <- c(Vasso, rep(1, ns_cox[[m]])) # derivative is 1 for linear time
                    }
                  }else{
                    if(assoCur %in% c("CV", "CV_CS")){
                      # evaluate f function of time at time points re.weight for current value association in survival
                      Vasso <- c(Vasso, unname(sapply(re.weight[[m]], paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == modelFE[[k]][[1]][j])))))
                    }
                    if(assoCur %in% c("CS", "CV_CS")){
                      # derivative of time function
                      Vasso <- c(Vasso, numDeriv::grad(get(paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == modelFE[[k]][[1]][j]))), re.weight[[m]]))
                    }
                  }
                }
                if(assoCur == "SRE"){
                  Vasso <- c(Vasso, rep(NA, ns_cox[[m]]))
                }
              }
              assoInfo <- rbind(assoInfo, c(assoCur, dim(data_cox[[m]])[1])) # keep record of what is done to avoid redoing it
              if(assoCur=="CV_CS") assoInfo <- rbind(assoInfo, c("CV", dim(data_cox[[m]])[1]), c("CS", dim(data_cox[[m]])[1]))
            }
          }
          assign(paste0(modelFE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), unname(modelFE[[k]][[2]][,j]), Vasso))
          # (rep NA NAvect is to fill the k-1 first markers where the marker k has no effect)
        }else{ # if association length is zero, there is no vector for association after the likelihood part for marker k
          assign(paste0(modelFE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), unname(modelFE[[k]][[2]][,j])))
        }
      }
      # add NA to match size of all markers until k for fixed effects of previous k-1 markers
      dataFE <- lapply(dataFE, function(x) append(x, rep(NA, dim(modelFE[[k]][[2]])[1]+length(Vasso))))
      tempNames <- names(dataFE) # save names before adding new items
      dataFE <- append(dataFE, mget(paste0(modelFE[[k]][[1]][1:length(modelFE[[k]][[1]])], "_L",k))) # add new items for marker k
      names(dataFE) <- c(tempNames, paste0(modelFE[[k]][[1]][1:length(modelFE[[k]][[1]])], "_L",k)) # set full vector of names

      for(j in 1:length(modelRE[[k]][[1]])){ # random effects
        Vasso <- NULL # vector for association part
        Wasso <- NULL # vector for association part (w = weight)
        if(length(assoc)!=0){
          assoInfoRE <- NULL
          assoCurRE <- NULL # current association
          for(m in 1:M){
            assoCurRE <- assoc[[k]][m]
            if(!ns_cox[[m]] %in% assoInfoRE[which(assoCurRE == assoInfoRE[,1]),2]){
              if(assoCurRE=="CV_CS"){ # if one of CV or CS is already done, don't redo it!
                if(ns_cox[[m]] %in% assoInfoRE[which("CV" == assoInfoRE[,1]),2]) assoCurRE <- "CS"
                if(ns_cox[[m]] %in% assoInfoRE[which("CS" == assoInfoRE[,1]),2]) assoCurRE <- "CV"
              }
              if(!(modelRE[[k]][[1]][j] %in% c(timeVar, c(paste0("f", 1:NFT, timeVar))))){ # if random effect j of marker k is not a time-dependent variable
                if(modelRE[[k]][[1]][j]=="Intercept"){
                  if(assoCurRE %in% c("CV", "CV_CS", "SRE")){
                    Vasso <- c(Vasso, c(IDre +  id_cox[[m]])) # individual id (unique for this vector of random effects)
                    Wasso <- c(Wasso, rep(1, ns_cox[[m]])) # weight is 1 because not time-dependent
                  }
                }else{
                  # establish id for a given variable
                  if(assoCurRE %in% c("CV", "CV_CS", "SRE", "SRE_ind")){
                    correspondID <- cbind(unique(data_cox[[m]][, modelRE[[k]][[1]][j]]), 1:length(unique(data_cox[[m]][, modelRE[[k]][[1]][j]])))
                    idVar <- unname(sapply(data_cox[[m]][, modelRE[[k]][[1]][j]], function(x) correspondID[which(correspondID[,1]==x),2])) #set id for random effect
                    if(assoCurRE != "SRE_ind"){
                      Vasso <- c(Vasso, c(IDre +  idVar)) # individual id (unique for this vector of random effects)
                      Wasso <- c(Wasso, rep(1, ns_cox[[m]])) # weight is 1 because not time-dependent
                    }
                  }
                }
                if(assoCurRE %in% c("CS", "CV_CS")){
                  Vasso <- c(Vasso, rep(NA, ns_cox[[m]]))
                  Wasso <- c(Wasso, rep(1, ns_cox[[m]]))
                }
              }else{
                if(modelRE[[k]][[1]][j] == timeVar){
                  if(assoCurRE %in% c("CV", "CV_CS", "SRE")){
                    Vasso <- c(Vasso, c(IDre +  id_cox[[m]])) # unique individual id
                    Wasso <- c(Wasso, re.weight[[m]]) # linear time weight
                  }
                  if(assoCurRE %in% c("CS", "CV_CS")){
                    Vasso <- c(Vasso, c(IDre +  id_cox[[m]])) # unique individual id
                    Wasso <- c(Wasso, rep(1, ns_cox[[m]])) # linear time weight
                  }
                }else{
                  # evaluate f function of time at time points re.weight for current value association in survival
                  if(assoCurRE %in% c("CV", "CV_CS")){
                    Vasso <- c(Vasso, c(IDre +  id_cox[[m]])) # unique individual id
                    Wasso <- c(Wasso, unname(sapply(re.weight[[m]], paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == modelRE[[k]][[1]][j])))))
                  }
                  if(assoCurRE %in% c("CS", "CV_CS")){
                    Vasso <- c(Vasso, c(IDre +  id_cox[[m]])) # unique individual id
                    Wasso <- c(Wasso, numDeriv::grad(get(paste0("f", which(c(paste0("f", 1:NFT, timeVar)) == modelRE[[k]][[1]][j]))), re.weight[[m]]))
                  }
                }
              }
              assoInfoRE <- rbind(assoInfoRE, c(assoCurRE, ns_cox[[m]]))
              if(assoCurRE=="CV_CS") assoInfoRE <- rbind(assoInfoRE, c("CV", ns_cox[[m]]), c("CS", ns_cox[[m]]))
            }
          }
        }
        if(!(modelRE[[k]][[1]][j] %in% c(timeVar, c(paste0("f", 1:NFT, timeVar))))){
          if(modelRE[[k]][[1]][j]=="Intercept"){
            assign(paste0("ID",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), IDre + as.integer(dataL[,id]), Vasso)) # assign variable with dynamic name for random effect
            assign(paste0("W",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), unname(modelRE[[k]][[2]][,j]), Wasso)) # assign variable with dynamic name for associated weight
          }else{
            if(exists("data_cox")){
              idVar <- unname(sapply(modelRE[[k]][[2]][,j], function(x) correspondID[which(correspondID[,1]==x),2])) #set id for random effect
              assign(paste0("ID",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), IDre + idVar, Vasso)) # assign variable with dynamic name for random effect
              assign(paste0("W",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), rep(1, length(idVar)), Wasso)) # assign variable with dynamic name for associated weight
            }else{
              correspondID <- cbind(unique(modelRE[[k]][[2]][,j]), 1:length(unique(modelRE[[k]][[2]][,j])))
              idVar <- unname(sapply(modelRE[[k]][[2]][,j], function(x) correspondID[which(correspondID[,1]==x),2])) #set id for random effect
            }
          }
        }else{
          assign(paste0("ID",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), IDre + as.integer(dataL[,id]), Vasso)) # assign variable with dynamic name for random effect
          assign(paste0("W",modelRE[[k]][[1]][j], "_L",k), c(rep(NA, NAvect), unname(modelRE[[k]][[2]][,j]), Wasso)) # assign variable with dynamic name for associated weight
        }
        if(!is.null(Vasso)){ # update the unique id counter so that it knows where to start at the next iteration
          IDre <- tail(na.omit(Vasso),1)
        } else{
          IDre <- tail(get(paste0("ID",modelRE[[k]][[1]][j], "_L",k)),1)
        }
      }
      assoRE <- NULL
      # set up outcome part (need to add association terms as outcomes too, equal to zero)
      outC <- list(modelYL[[k]][[2]]) # outcome values for marker k
      names(outC) <- modelYL[[k]][[1]] # name
      fam <- c(fam, family[[k]])
      famCtrl <- append(famCtrl, list(list(link=link[k])))
      if(length(assoc)!=0){ # set up the association part
        uv <- c(rep(NA, length(dataL[,id])+NAvect)) # used if current value association
        wv <- c(rep(NA, length(dataL[,id])+NAvect)) # corresponding weight
        us <- c(rep(NA, length(dataL[,id])+NAvect)) # used if current slope association
        ws <- c(rep(NA, length(dataL[,id])+NAvect)) # corresponding weight
        usre <- c(rep(NA, length(dataL[,id])+NAvect)) # used if shared random effect association
        wsre <- c(rep(NA, length(dataL[,id])+NAvect)) # corresponding weight
        assoInfo2 <- NULL
        assoCur2 <- NULL
        assoC <- list() # outcome (association)
        assoInfo3 <- NULL
        NAasso <- 0 # NA in case of multiple associations
        for(m in 1:M){ # for each unique association term for marker k
          assoCur2 <- assoc[[k]][m]
          if(!dim(data_cox[[m]])[1] %in% assoInfo2[which(assoCur2 == assoInfo2[,1]),2]){
            if(assoCur2=="CV_CS"){ # if one of CV or CS is already done, don't redo it!
              if(ns_cox[[m]] %in% assoInfo2[which("CV" == assoInfo2[,1]),2]) assoCur2 <- "CS"
              if(ns_cox[[m]] %in% assoInfo2[which("CS" == assoInfo2[,1]),2]) assoCur2 <- "CV"
            }
            if(assoCur2 %in% c("CV", "CV_CS")){ # current value
              uv <- c(uv, unlist(data_cox[[m]][paste0("CV_L", k, "_S", m)])) # get the id corresponding to this association to match it
              wv <- c(wv, rep(-1, ns_cox[[m]]))
              us <- c(us, rep(NA, ns_cox[[m]]))
              ws <- c(ws, rep(NA, ns_cox[[m]]))
              usre <- c(usre, rep(NA, ns_cox[[m]]))
              wsre <- c(wsre, rep(NA, ns_cox[[m]]))
              # fill association part with NA for the part of the vector corresponding to the likelihood of the marker (where we set the outcome values)
              outC[[1]] <- c(outC[[1]], rep(NA, ns_cox[[m]])) # outcome part is NA
              assoC <- append(assoC, list(c(rep(NA, length(dataL[, id])+NAvect+NAasso), rep(0, ns_cox[[m]])))) # set NA for the marker's likelihood part
              names(assoC)[length(assoC)] <- paste0("CV_L", k, "_S", m, m) # add dynamic name to make it unique
              assoInfo3 <- rbind(assoInfo3, c("CV", ns_cox[[m]]))
              if(length(assoC)>1){
                for(tassoc in 1:(length(assoC)-1)){
                  assoC[[tassoc]] <- c(assoC[[tassoc]], rep(NA, ns_cox[[m]]))
                }
              }
              NAasso <- NAasso + ns_cox[[m]]
              fam <- c(fam, "gaussian")
              famCtrl <- append(famCtrl, list(list(hyper = list(prec = list(initial = 12, fixed=TRUE)))))
            }
            if(assoCur2 %in% c("CS", "CV_CS")){ # current slope
              us <- c(us, unlist(data_cox[[m]][paste0("CS_L", k, "_S", m)]))
              ws <- c(ws, rep(-1, ns_cox[[m]]))
              uv <- c(uv, rep(NA, ns_cox[[m]]))
              wv <- c(wv, rep(NA, ns_cox[[m]]))
              usre <- c(usre, rep(NA, ns_cox[[m]]))
              wsre <- c(wsre, rep(NA, ns_cox[[m]]))
              outC[[1]] <- c(outC[[1]], rep(NA, ns_cox[[m]]))
              assoC <- append(assoC, list(c(rep(NA, length(dataL[, id])+NAvect+NAasso), rep(0, ns_cox[[m]]))))
              names(assoC)[length(assoC)] <- paste0("CS_L", k, "_S", m, m)
              assoInfo3 <- rbind(assoInfo3, c("CS", ns_cox[[m]]))
              if(length(assoC)>1){
                for(tassoc in 1:(length(assoC)-1)){
                  assoC[[tassoc]] <- c(assoC[[tassoc]], rep(NA, ns_cox[[m]]))
                }
              }
              NAasso <- NAasso + ns_cox[[m]]
              fam <- c(fam, "gaussian")
              famCtrl <- append(famCtrl, list(list(hyper = list(prec = list(initial = 12, fixed=TRUE)))))
            }
            if(assoCur2 %in% c("SRE")){ # individual deviation
              usre <- c(usre, unlist(data_cox[[m]][paste0("SRE_L", k, "_S", m)]))
              wsre <- c(wsre, rep(-1, ns_cox[[m]]))
              uv <- c(uv, rep(NA, ns_cox[[m]]))
              wv <- c(wv, rep(NA, ns_cox[[m]]))
              us <- c(us, rep(NA, ns_cox[[m]]))
              ws <- c(ws, rep(NA, ns_cox[[m]]))
              outC[[1]] <- c(outC[[1]], rep(NA, ns_cox[[m]]))
              assoC <- append(assoC, list(c(rep(NA, length(dataL[, id])+NAvect+NAasso), rep(0, ns_cox[[m]]))))
              names(assoC)[length(assoC)] <- paste0("SRE_L", k, "_S", m, m)
              assoInfo3 <- rbind(assoInfo3, c("SRE", ns_cox[[m]]))
              if(length(assoC)>1){
                for(tassoc in 1:(length(assoC)-1)){
                  assoC[[tassoc]] <- c(assoC[[tassoc]], rep(NA, ns_cox[[m]]))
                }
              }
              NAasso <- NAasso + ns_cox[[m]]
              fam <- c(fam, "gaussian")
              famCtrl <- append(famCtrl, list(list(hyper = list(prec = list(initial = 12, fixed=TRUE)))))
            }
            assoInfo2 <- rbind(assoInfo2, c(assoCur2, ns_cox[[m]]))
            if(assoCur2=="CV_CS") assoInfo2 <- rbind(assoInfo2, c("CV", ns_cox[[m]]), c("CS", ns_cox[[m]]))
          }
        }
        YL <- lapply(YL, function(x) append(x, rep(NA, length(outC[[1]])))) # add NA to match size of all markers until k
        outC[[1]] <- c(rep(NA, NAvect), outC[[1]])
        YL <- append(YL, c(outC, assoC)) # add outcome and association

        assign(paste0("uv",k), unname(uv)) # assign association with dynamic variable name
        assign(paste0("wv",k), wv) # associated weight
        if("CV" %in% assoc[[k]] | "CV_CS" %in% assoc[[k]]) assoRE <- c(assoRE, paste0("uv",k), paste0("wv",k)) # random effects association
        assign(paste0("us",k), unname(us))
        assign(paste0("ws",k), ws)
        if("CS" %in% assoc[[k]] | "CV_CS" %in% assoc[[k]]) assoRE <- c(assoRE, paste0("us",k), paste0("ws",k))
        assign(paste0("usre",k), unname(usre))
        assign(paste0("wsre",k), wsre)
        if("SRE" %in% assoc[[k]]) assoRE <- c(assoRE, paste0("usre",k), paste0("wsre",k))
      }else{ # if no association, directly add the marker's likelihood part without having to set up an association part
        YL <- lapply(YL, function(x) append(x, rep(NA, length(outC[[1]]))))
        outC[[1]] <- c(rep(NA, NAvect), outC[[1]])
        YL <- append(YL, outC) # add outcome and association
      }
      NAvect <- length(YL[[k]]) # size of the vector of NA for next iteration

      # add NA to match size of all markers until k for fixed effects of previous k-1 markers
      dataRE <- lapply(dataRE, function(x) append(x, rep(NA, dim(modelRE[[k]][[2]])[1]+length(Vasso))))
      tempNames <- names(dataRE)
      dataRE <- append(dataRE, mget(c(paste0("ID",modelRE[[k]][[1]][1:length(modelRE[[k]][[1]])], "_L",k),
                                      paste0("W",modelRE[[k]][[1]][1:length(modelRE[[k]][[1]])], "_L",k),
                                      assoRE))) # random effects data part (with likelihoof for marker k followed by association for marker k)
      names(dataRE) <- c(tempNames, paste0("ID",modelRE[[k]][[1]][1:length(modelRE[[k]][[1]])], "_L",k),
                         paste0("W",modelRE[[k]][[1]][1:length(modelRE[[k]][[1]])], "_L",k),
                         assoRE)
      dataRE <- lapply(dataRE, function(x) unname(x)) # clean data: remove useless names of some parts of the vectors
    }

    REstruc=NULL # store random effects structure for summary()
    REstruc1 <- sapply(modelRE,"[[",1)
    if(class(REstruc1)[1]=="matrix"){
      for(k in 1:ncol(REstruc1)){
        for(l in 1:length(REstruc1[, k])){
          REstruc <- c(REstruc, paste0(REstruc1[, k][l], "_L", k))
        }
      }
    }else{
      for(k in 1:length(REstruc1)){
        for(l in 1:length(REstruc1[[k]])){
          REstruc <- c(REstruc, paste0(REstruc1[[k]][l], "_L", k))
        }
      }
    }
  }
  ################################################################## joint fit
  if(is_Long) jointdf = data.frame(dataFE, dataRE, YL) # dataset with fixed and random effects as well as outcomes for the K markers
  # at this stage all the variables have unique mname that refers to the number of the marker (k) ot the number of the survival outcome (m)
  if(is_Surv){
    if(is_Long){
      joint.data <- c(as.list(inla.rbind.data.frames(jointdf, Map(c,data_cox[1:M]))), dlCox)
      Yjoint = as.list(rbind.data.frame(joint.data[c(names(YL))]))
      for(m in 1:M){
        if(basRisk[[m]] %in% c("rw1", "rw2")){
          Yjoint <- append(Yjoint, joint.data[paste0("y", m, "..coxph")])
        }else{
          if(!is.null(assoc)) YS_assoc <- unlist(assoc[1:K])[seq(m, K*M, by=M)] else YS_assoc <- NULL # extract K association terms associated to time-to-event m
          if(!TRUE %in% c(c("CV", "CS", "CV_CS", "SRE") %in% YS_assoc)){
            assign(paste0(get(paste0("cox_event_", m))$formula[[2]]), list(inla.surv(time = joint.data[[paste0("surv", m, "time")]], event = joint.data[[paste0("y", m, "..coxph")]])))
            Yjoint <- append(Yjoint, get(get(paste0("cox_event_", m))$formula[[2]]))
          }else{
            assign(paste0(get(paste0("cox_event_", m))$formula[[2]]), inla.surv(truncation = joint.data[[paste0("baseline", m, ".hazard.time")]],
                                                                                time = joint.data[[paste0("baseline", m, ".hazard.time")]] + joint.data[[paste0("E..coxph")]],
                                                                                event = joint.data[[paste0("y", m, "..coxph")]]))
            Yjoint <- append(Yjoint, list(get(get(paste0("cox_event_", m))$formula[[2]])))
          }
        }
      }
      joint.data$Yjoint <- Yjoint # all the data is in this object (longitudinal, survival)
    }else{
      joint.data <- c(as.list(inla.rbind.data.frames(Map(c,data_cox[1:M]))), dlCox)
      Yjoint <- NULL
      for(m in 1:M){
        if(basRisk[[m]] %in% c("rw1", "rw2")){
          Yjoint <- append(Yjoint, joint.data[paste0("y", m, "..coxph")])
        }else{
          if(!is.null(assoc)) YS_assoc <- unlist(assoc[1:K])[seq(m, K*M, by=M)] else YS_assoc <- NULL # extract K association terms associated to time-to-event m
          if(!TRUE %in% c(c("CV", "CS", "CV_CS", "SRE") %in% YS_assoc)){
            assign(paste0(get(paste0("cox_event_", m))$formula[[2]]), list(inla.surv(time = joint.data[[paste0("surv", m, "time")]], event = joint.data[[paste0("y", m, "..coxph")]], cure=get(paste0("cure_",m)))))
            Yjoint <- append(Yjoint, get(get(paste0("cox_event_", m))$formula[[2]]))
          }else{
            assign(paste0(get(paste0("cox_event_", m))$formula[[2]]), inla.surv(time = joint.data[[paste0("baseline", m, ".hazard.time")]],
                                                                                time2 = joint.data[[paste0("baseline", m, ".hazard.time")]] + joint.data[[paste0("E..coxph")]],
                                                                                event = joint.data[[paste0("y", m, "..coxph")]],
                                                                                cure=get(paste0("cure_",m))))
            Yjoint <- append(Yjoint, list(get(get(paste0("cox_event_", m))$formula[[2]])))
          }
        }
      }
      joint.data$Yjoint <- Yjoint # all the data is in this object (survival)
    }
    # formula: survival part
    if(M>1){
      for(m in 1:M){ # if more than one survival submodel then we need to merge the formulas
        if(m!=M){
          formAdd <- paste0("Yjoint ~ . + ", strsplit(as.character(get(paste0("formS", m+1))), "~")[[3]])
          if(m==1) FormAct <- get(paste0("formS", m)) else FormAct <- formulaSurv
          formulaSurv = update(FormAct, formAdd) # update to have a unique formula for survival outcomes up to m
        }
        if(!is.null(modelYS[[m]]$RE_matS)){ # random effects in survival model m
          formulaSurv = update(formulaSurv, formAddS[[m]])
        }
      }
    }else{
      if(!is.null(modelYS[[m]]$RE_matS)){ # random effects in survival model m
        formulaSurv = update(get(paste0("formS", 1:M)), formAddS[[m]])
      }else{
        formulaSurv = update(get(paste0("formS", 1:M)), "Yjoint ~ .") # if only one survival outcome, directly extract the corresponding formula
      }
    }
  }else{
    joint.data <- as.list(jointdf) # remove Y not used here?
    Yjoint = as.list(rbind.data.frame(joint.data[names(YL)]))
    joint.data$Yjoint <- Yjoint # all the data is in this object (longitudinal, survival)
  }

  if(is_Long){
    # formula: random effects part
    formulaRand <- vector("list", K) # model for longitudinal markers
    if(corLong){
      nTot <- 0
      sTot <- 0
      for(k in 1:K){
        nTot <- nTot + length(modelRE[[k]][[1]]) # get number of random effects if they are correlated
        sTot <- sTot + max(unlist(Nid)) * length(modelRE[[k]][[1]]) # get the sum of the sizes # Nid[[k]]
      }
      if(nTot>10) stop(paste0("The maximum number of correlated random effects is 10 and you request ", nTot,
                              ". Please reduce the number of random effects or assume independent longitudinal
                              markers (set parameter corLong to FALSE). If you need to overcome this limit,
                              please contact us (INLAjoint@gmail.com)."))
    }

    # formula: association part
    formulaAssoc <- vector("list", K) # model for longitudinal markers
    for(k in 1:K){ # for each marker k
      form1 <- NULL
      form2 <- NULL
      if(corLong){
        if(K==1) stop("Only one longitudinal marker is detected but 'corLong' is set to TRUE. Please set 'corLong' to FALSE
                       or include at least two longitudinal markers to have their random effects correlated.")
        if(k==1){
          if(length(modelRE[[k]][[1]])==1){
            form1 <- paste("f(", paste0("ID",modelRE[[k]][[1]], "_L",k)[1],",", paste0("W",modelRE[[k]][[1]], "_L",k)[1],", model = 'iidkd', order=",nTot,
                           ", n =", sTot,", constr = F, hyper = list(theta1 = list(param = c(", control$PriorRandom$r,", ", paste(c(rep(control$PriorRandom$R, nTot), rep(0, (nTot*nTot-nTot)/2)), collapse=","), "))))")
          }else if(length(modelRE[[k]][[1]])>1){ # if two random effects, use cholesky parameterization (i.e., iidkd)
            form1 <- paste("f(", paste0("ID",modelRE[[k]][[1]], "_L",k)[1],",", paste0("W",modelRE[[k]][[1]], "_L",k)[1],", model = 'iidkd',
                order = ",nTot,", n =", sTot,", constr = F, hyper = list(theta1 = list(param = c(", control$PriorRandom$r,", ", paste(c(rep(control$PriorRandom$R, nTot), rep(0, (nTot*nTot-nTot)/2)), collapse=","), "))))")
            for(fc in 2:length(modelRE[[k]][[1]])){
              form2 <- paste(c(form2, paste0("f(",paste0("ID",modelRE[[k]][[1]], "_L",k)[fc],",", paste0("W",modelRE[[k]][[1]], "_L",k)[fc],",
                                  copy = ",paste0("'ID",modelRE[[1]][[1]], "_L1'")[1],")")), collapse="+")
            }
          }
        }else{
          for(fc in 1:length(modelRE[[k]][[1]])){
            form2 <- paste(c(form2, paste0("f(",paste0("ID",modelRE[[k]][[1]], "_L",k)[fc],",", paste0("W",modelRE[[k]][[1]], "_L",k)[fc],",
                                  copy = ",paste0("'ID",modelRE[[1]][[1]], "_L1'")[1],")")), collapse="+")
          }
        }

      }else{
        if(length(modelRE[[k]][[1]])==1){ # if only one random effect, need to use "iid"
          form1 <- paste("f(", paste0("ID",modelRE[[k]][[1]], "_L",k)[1],",", paste0("W",modelRE[[k]][[1]], "_L",k)[1],", model = 'iid',
                n =", Nid[[k]] * length(modelRE[[k]][[1]]),", constr = F)")
        }else if(length(modelRE[[k]][[1]])>1){ # if two random effects, use cholesky parameterization (i.e., iidkd)
          form1 <- paste("f(", paste0("ID",modelRE[[k]][[1]], "_L",k)[1],",", paste0("W",modelRE[[k]][[1]], "_L",k)[1],", model = 'iidkd',
                 order = ",length(modelRE[[k]][[1]]),", n =", Nid[[k]] * length(modelRE[[k]][[1]]),", constr = F, hyper = list(theta1 =
                 list(param = c(", control$PriorRandom$r,", ", paste(c(rep(control$PriorRandom$R, length(modelRE[[k]][[1]])), rep(0, (length(modelRE[[k]][[1]])*length(modelRE[[k]][[1]])-length(modelRE[[k]][[1]]))/2)), collapse=","), "))))")
          for(fc in 2:length(modelRE[[k]][[1]])){
            form2 <- paste(c(form2, paste0("f(",paste0("ID",modelRE[[k]][[1]], "_L",k)[fc],",", paste0("W",modelRE[[k]][[1]], "_L",k)[fc],",
                                  copy = ",paste0("'ID",modelRE[[k]][[1]], "_L",k, "'")[1],")")), collapse="+")
          }
        }
      }
      formulaRand[[k]] <- paste(c(form1, form2), collapse="+")

      if(length(assoc)!=0){# if there is at least one association term
        for(Nassoc in 1:length(assoc[[k]])){ # for each association term included for marker k (should have length M)
          if("CV" == assoc[[k]][[Nassoc]]){ # if current value
            if(!is.null(formulaAssoc[[k]])){ # if there is something in this object, first verify that the current value is not already set up for this marker
              if(gsub("\\s", "", paste0("f(uv",k,", wv",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]){
                # if current value is already set up for this marker, just add the association part
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")
              }else{
                # if current value is not set up for this marker, first set it up (with uv and wv) and then add the association part (CV_L)
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(uv",k,", wv",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                         paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }
            }else{

              formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(uv",k,", wv",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                       paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
            }
          }else if("CS" == assoc[[k]][[Nassoc]]){ # current slope
            if(!is.null(formulaAssoc[[k]])){
              if(gsub("\\s", "", paste0("f(us",k,", ws",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")
              }else{

                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(us",k,", ws",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),

                                                                         paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }
            }else{
              formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(us",k,", ws",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                       paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
            }
          }else if("SRE" == assoc[[k]][[Nassoc]]){ # shared random effects
            if(!is.null(formulaAssoc[[k]])){
              if(gsub("\\s", "", paste0("f(usre",k,", wsre",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("usre",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")
              }else{

                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(usre",k,", wsre",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),

                                                                         paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("usre",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }
            }else{
              formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(usre",k,", wsre",k,", model = 'iid',  hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                       paste0("f(", paste0(assoc[[k]][Nassoc], "_L", k, "_S", Nassoc), ", copy='", paste0("usre",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
            }
          }else if("SRE_ind" == assoc[[k]][[Nassoc]]){ # shared random effects independent
            for(i in 1:length(modelRE[[k]][[1]])){
              formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste0("f(SRE_",modelRE[[k]][[1]][i] , "_L", k, "_S", Nassoc, ", copy='", paste0("ID",paste0(modelRE[[k]][[1]][i]),"_L", k,"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"))), collapse="+")
            }
          }else if("CV_CS" == assoc[[k]][[Nassoc]]){ # current value + current slope
            if(!is.null(formulaAssoc[[k]])){
              if(!(gsub("\\s", "", paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]] & paste0("f(us",k,", ws",k,", hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)") %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]])){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                         paste0("f(", paste0("CV_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"),
                                                                         paste0("f(us",k,", ws",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                         paste0("f(", paste0("CS_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")

              }else if(!(gsub("\\s", "", paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]) & paste0("f(us",k,", ws",k,", hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)") %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                         paste0("f(", paste0("CV_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"),
                                                                         paste0("f(", paste0("CS_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }else if(gsub("\\s", "", paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]] & !(paste0("f(us",k,", ws",k,", hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)") %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]])){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(", paste0("CV_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"),
                                                                         paste0("f(us",k,", ws",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                         paste0("f(", paste0("CS_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }else if(gsub("\\s", "", paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)")) %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]] & paste0("f(us",k,", ws",k,", hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)") %in% strsplit(gsub("\\s", "", formulaAssoc[[k]]), "\\+")[[1]]){
                formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(", paste0("CV_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"),
                                                                         paste0("f(", paste0("CS_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
              }
            }else{
              formulaAssoc[[k]] <-  paste(c(formulaAssoc[[k]], paste(c(paste0("f(uv",k,", wv",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                       paste0("f(", paste0("CV_L", k, "_S", Nassoc), ", copy='", paste0("uv",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))"),



                                                                       paste0("f(us",k,", ws",k,", model = 'iid', hyper = list(prec = list(initial = -6,fixed = TRUE)), constr = F)"),
                                                                       paste0("f(", paste0("CS_L", k, "_S", Nassoc), ", copy='", paste0("us",k),"', hyper = list(beta = list(fixed = FALSE,param = c(", control$priorAssoc$mean,",", control$priorAssoc$prec,"), initial = ", assocInit, ")))")), collapse="+")), collapse="+")
            }
          }
        }
      }else formulaAssoc <- NULL
    }
    # formula: longitudinal part
    # merge outcome, fixed effects, random effects and association terms
    formulaLong = formula(paste(c("Yjoint ~ . -1", names(dataFE), formulaRand[1:K], formulaAssoc[1:K]), collapse="+"))
  }else{
    formulaAssoc <- NULL
    formulaLong <- NULL
  }

  # final formula
  # merge survival and longitudinal parts formulas
  if(is_Long & is_Surv){
    formulaJ <- update(formulaSurv, formulaLong)
  }else if(is_Long & !is_Surv){
    formulaJ <- formulaLong
  }else if(!is_Long & is_Surv){
    formulaJ <- formulaSurv
  }
  if(is_Long){
    for(k in 1:K){
      if("poisson" == family[[k]]){ # if longitudinal marker k is poisson, need to set up the E equal to 1 for the part of the vector corresponding to this marker
        if(is_Surv){
          joint.data$E..coxph[which(!is.na(joint.data$Yjoint[[which(names(joint.data$Yjoint) == modelYL[[k]][[1]])]]))] <- 1
        }else{
          joint.data <- append(joint.data, list(c(ifelse(is.na(joint.data$Yjoint[[which(names(joint.data$Yjoint) == modelYL[[k]][[1]])]]), NA, 1))))
          names(joint.data)[length(names(joint.data))] <- "E..coxph"
        }
      }else if("binomial" == family[[k]]){ # for binomial, make sure we have integers)
        joint.data$Yjoint[[which(names(joint.data$Yjoint) == modelYL[[k]][[1]])]] <- as.integer(as.factor(joint.data$Yjoint[[which(names(joint.data$Yjoint) == modelYL[[k]][[1]])]]))-1
#        linkBinom <- which(names(joint.data$Yjoint) == modelYL[[k]][[1]])
      }
    }
    if(length(assoc)!=0){ # for the association terms, we have to add the gaussian family and specific hyperparameters specifications
      if(is_Surv){
        familySurv <- NULL
        for(m in 1:M){
          familySurv <- c(familySurv, ifelse(basRisk[[m]] %in% c("rw1", "rw2"), "poisson", basRisk[[m]]))
          if(!is.null(get(paste0("cure_",m)))){ # set priors for all components of the cure model
            HYP <- "list("
            for(i in 1:dim(get(paste0("cure_",m)))[2]){
              HYP <- paste0(HYP, "beta",i,"= list(prior = 'normal', param = c(", control$priorFixed$mean, ", ", control$priorFixed$prec, "))")
              if(i != dim(get(paste0("cure_",m)))[2]) HYP <- paste0(HYP, ", ")
            }
            HYP <- paste0(HYP, ")")
            if(i==dim(get(paste0("cure_",m)))[2]) HYPER <- eval(parse(text=HYP))
          } else HYPER <- list()
          HYPER <- append(HYPER, control$baselineHyper)
          famCtrl <- append(famCtrl, ifelse(basRisk[[m]]=="weibullsurv", list(list(variant=variant, hyper=HYPER)), list(list())))
        }
        fam <- unlist(c(fam, familySurv))
      }
    }else if(is_Surv){ # if no association directly add the survival part
      familySurv <- NULL
      famCtrl <- NULL
      for(m in 1:M){
        familySurv <- c(familySurv, ifelse(basRisk[[m]] %in% c("rw1", "rw2"), "poisson", basRisk[[m]]))
        if(!is.null(get(paste0("cure_",m)))){ # set priors for all components of the cure model
          HYP <- "list("
          for(i in 1:dim(get(paste0("cure_",m)))[2]){
            HYP <- paste0(HYP, "beta",i,"= list(prior = 'normal', param = c(", control$priorFixed$mean, ", ", control$priorFixed$prec, "))")
            if(i != dim(get(paste0("cure_",m)))[2]) HYP <- paste0(HYP, ", ")
          }
          HYP <- paste0(HYP, ")")
          if(i==dim(get(paste0("cure_",m)))[2]) HYPER <- eval(parse(text=HYP))
        } else HYPER <- list()
        HYPER <- append(HYPER, control$baselineHyper)
        famCtrl <- append(famCtrl, ifelse(basRisk[[m]]=="weibullsurv", list(list(variant=variant, hyper=HYPER)), list(list())))
      }
      fam <- unlist(c(family, familySurv))
    }else if(!is_Surv){
      # fam <- unlist(c(family))
      # for(k in 1:K){
      #   famCtrl <- c(famCtrl, list(list(link=link[k])))
      # }
    }
  }else{
    fam <- NULL
    famCtrl <- NULL
    for(m in 1:M){
      fam <- c(fam, ifelse(basRisk[[m]] %in% c("rw1", "rw2"), "poisson", basRisk[[m]]))
      if(!is.null(get(paste0("cure_",m)))){ # set priors for all components of the cure model
        HYP <- "list("
        for(i in 1:dim(get(paste0("cure_",m)))[2]){
          HYP <- paste0(HYP, "beta",i,"= list(prior = 'normal', param = c(", control$priorFixed$mean, ", ", control$priorFixed$prec, "))")
          if(i != dim(get(paste0("cure_",m)))[2]) HYP <- paste0(HYP, ", ")
        }
        HYP <- paste0(HYP, ")")
        if(i==dim(get(paste0("cure_",m)))[2]) HYPER <- eval(parse(text=HYP))
      } else HYPER <- list()
      HYPER <- append(HYPER, control$baselineHyper)
      famCtrl <- append(famCtrl, ifelse(basRisk[[m]]=="weibullsurv", list(list(variant=variant, hyper=HYPER)), list(list())))
    }
  }
  RMVN <- NULL # "ReMoVe Names" : for random walks, we remove the intercept and the unconstrained random walk will give it
  if(is_Surv){
    for(m in 1:M){
      if(basRisk[m]%in%c("rw1", "rw2")){
        RMVN <- c(RMVN, paste0("Intercept_S", m))
      }
    }
  }
  if(is.null(control$cmin)) control$cmin <- ifelse(is.null(unlist(cureVar)), 0,-Inf) # for inla(), needs to be set to 1 if a mixture cure component is in the model for stability
  # if no survival component, need to remove dot in formula
  if(!is_Surv) formulaJ <- formula(paste("Yjoint~-1", strsplit(as.character(formulaJ)[3], "\\. - 1")[[1]][2]))
  # fix issue with formula
  formulaJ <- eval(parse(text=paste0(as.character(formulaJ)[2], as.character(formulaJ)[1], as.character(formulaJ)[3])))
  if(length(joint.data$Yjoint)==1) joint.data$Yjoint <- joint.data$Yjoint[[1]]
  res <- inla(formulaJ,family = fam,
              data=joint.data,
              control.fixed = list(mean=control$priorFixed$mean, prec=control$priorFixed$prec,
                                   mean.intercept=control$priorFixed$mean.intercept, prec.intercept=control$priorFixed$prec.intercept, remove.names=RMVN),
              control.family = famCtrl, inla.mode = "experimental",
              control.compute=list(config = cfg, dic=T, waic=T, cpo=cpo,
                                   control.gcpo = list(enable = cpo,
                                                       num.level.sets = -1,
                                                       correct.hyperpar = TRUE)),
              E = joint.data$E..coxph,
              control.inla = list(int.strategy=int.strategy, cmin=control$cmin),#parallel.linesearch=T, cmin = 0
              safe=safemode, verbose=verbose, keep = keep)
  CLEANoutput <- c('summary.lincomb','mfarginals.lincomb','size.lincomb',
                   'summary.lincomb.derived','marginals.lincomb.derived','size.lincomb.derived','offset.linear.predictor',
                   'model.spde2.blc','summary.spde2.blc','marginals.spde2.blc','size.spde2.blc','model.spde3.blc','summary.spde3.blc',
                   'marginals.spde3.blc','size.spde3.blc','logfile','Q','graph','ok','model.matrix')
  res[CLEANoutput] <- NULL
  if(is_Surv) res$cureVar <- cureVar
  if(is_Long) res$famLongi <- unlist(family)
  if(exists("REstruc")) res$REstruc <- REstruc
  if(exists("REstrucS")) res$REstrucS <- REstrucS
  res$basRisk <- basRisk
  class(res) <- c("INLAjoint", "inla")
  return(res)
}


