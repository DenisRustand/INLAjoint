#' Setup survival part for outcome m
#' @description Setup survival part for outcome m (internal function)
#' input:
#' @param formula with inla.surv() object as an outcome
#' @param formLong formula from the longitudinal part, if any
#' @param dataSurv dataset(s) for the survival part
#' @param LSurvdat dataset for the longitudinal part converted to survival format (internal, used to get covariates
#' if missing in the survival dataset when sharing linear predictors including covariates from longitudinal into survival)
#' @param timeVar names of the variables that are time-dependent (only linear for now)
#' @param assoc association parameters between K longitudinal outcomes and M survival outcomes (list of K vectors of size M)
#' @param id name of the variable that gives the individual id
#' @param m identifies the outcome among 1:M time-to-event outcomes
#' @param K number of longitudinal outcomes
#' @param M number of survival outcomes
#' @param NFT maximum number of functions of time (fixed value)
#' @param corLong boolean that indicates if random effects across longitudinal markers are correlated,
#' when multiple longitudinal markers are included in the model
#' @param dataOnly boolean for internal use, indicates if only preparing data (i.e., not fitting the model)
#' @param SurvInfo information about survival submodels for internal use
#' @param strata covariate for stratified proportional hazards model
#' output:
#' @return YS_data includes the values of the survival outcome and covariates associated to this survival part,
#'         with the association parameters but the provided id are temporary and they will be updated after
#'         the cox expansion to make them unique and allow for time dependency
#' @return YSformF formula for this survival outcome (not including association parameters)
#' @importFrom lme4 nobars findbars
#' @export


setup_S_model <- function(formula, formLong, dataSurv, LSurvdat, timeVar, assoc, id, m, K, M, NFT, corLong, dataOnly, SurvInfo, strata){
  # Event outcome
  YS <- strsplit(as.character(formula), split="~")[[2]]
  FE_formS <- nobars(formula)
  RES <- findbars(formula) # random effects included? (i.e., frailty)
  # YS_FE_elements <- gsub("\\s", "", strsplit(strsplit(as.character(FE_formS), split="~")[[3]], split=c("\\+"))[[1]])
  FML <- ifelse(strsplit(as.character(FE_formS), split="~")[[3]]=="-1", 1, strsplit(as.character(FE_formS), split="~")[[3]])
  YSform2 <- formula(paste(" ~ ", FML))
  # if(length(which(sapply(dataSurv, class)=="factor"))>0){ # deal with factors when modalities are missing
  #   factors_columns <- which(sapply(dataSurv, class)=="factor")
  #   for(fctc in factors_columns){
  #     if(length(unique(dataSurv[, fctc]))< length(levels(dataSurv[, fctc]))) dataSurv[, fctc] <- as.integer(dataSurv[, fctc])-1
  #   }
  # }
  # YS_FE_elements2 <- YS_FE_elements[YS_FE_elements!="1" & YS_FE_elements!="-1"]
  if(length(all.vars(YSform2))>0){
    sapply(all.vars(YSform2), function(x) if(!(x %in% colnames(dataSurv))) stop(paste0("Covariate `", x, "` not found in survival dataset!")))
  }
  DFS <- model.matrix(YSform2, model.frame(YSform2, dataSurv, na.action=na.pass))
  if(colnames(DFS)[1]=="(Intercept)") colnames(DFS)[1] <- "Intercept"
  if(grepl("inla.surv", YS)){
    # attach(dataSurv)
    if(dataOnly & length(paste0(SurvInfo$survOutcome))>1) assign(paste0(SurvInfo$survOutcome)[2], dataSurv)
    YS <- with(dataSurv, eval(parse(text=YS)))
    YSname <- paste0("S", m)
    # detach(dataSurv)
  }else{
    YSname <- YS
    YS <- get(YS)
  }
  YS_data <- c(list(YS), as.list(as.data.frame(DFS)))
  names(YS_data)[1] <- YSname
  names(YS_data) <- paste0(gsub(":", ".X.", gsub("\\s", ".", names(YS_data))), "_S", m)
  YSformF <- formula(paste0(YSname, "_S", m, " ~ -1 +", paste0(paste0(gsub(":", ".X.", gsub("\\s", ".", colnames(DFS))), "_S", m, collapse="+"))))
  CLid <- 0 # keep track for unique id if corLong is true for shared random effects independently
  # strata?
  if(!is.null(strata)){
    YS_data <- append(YS_data, list(dataSurv[, strata]))
    names(YS_data)[length(YS_data)] <- strata
  }
  if(dataOnly & !is.null(id)){
    if(all.equal(dataSurv[[id]], LSurvdat[[id]])!=TRUE) LSurvdat[[id]] <- dataSurv[[id]]
  }
  # association
  if(length(assoc)!=0){
    YS_assoc <- unlist(assoc[1:K])[seq(m, K*M, by=M)] # extract K association terms associated to time-to-event m
    for(k in 1:length(YS_assoc)){
      if(TRUE %in% (YS_assoc %in% c("CV", "CS", "CS2", "CV_CS"))){
        # add covariates that are being shared through the association
        FE_form <- nobars(formLong[[k]])
        if(dim(LSurvdat)[1] != dim(dataSurv)[1] & !F %in% c(rownames(attr(terms(FE_form), "factors")) %in% colnames(dataSurv))){
          DFS2 <- as.data.frame(model.matrix(FE_form, model.frame(FE_form, dataSurv, na.action=na.pass)))
        }else{
          DFS2 <- as.data.frame(model.matrix(FE_form, model.frame(FE_form, LSurvdat[which(LSurvdat[[id]] %in% dataSurv[[id]]),], na.action=na.pass)))
        }
        removeVar <- NULL
        for(rmtvar in 1:length(strsplit(colnames(DFS2), ":"))){ # remove any component that contains a timeVar because it will not be useful
          if(TRUE %in% (c(timeVar, c(paste0("f", 1:NFT, "(", timeVar, ")"))) %in% strsplit(colnames(DFS2), ":")[[rmtvar]]) |
             TRUE %in% (names(YS_data) %in% strsplit(colnames(DFS2), ":")[[rmtvar]])) removeVar <- c(removeVar, rmtvar)
        }
        colNvar <- colnames(DFS2)
        DFS2 <- DFS2[,-c(which(colNvar %in% c("(Intercept)")), removeVar)]
        if(is.null(dim(DFS2))){
          YS_data <- append(YS_data, list(DFS2))
          names(YS_data)[length(names(YS_data))] <- colNvar[-c(which(colNvar %in% c("(Intercept)")), removeVar)]
        }else if(dim(DFS2)[2]>0){
          names(DFS2) <- gsub(":", ".X.", gsub("\\s", ".", names(DFS2)))
          YS_data <- append(YS_data, DFS2)
        }
        RW <- which(sapply(YS_data, length)[-1] != dim(dataSurv)[1])
        # overwrite variables that are not properly captured
        if(length(RW)>0 & (FALSE %in% c(YS_assoc=="CS"))){
          for(RWi in 1:length(RW)){
            if(length(unlist(sapply(colnames(dataSurv), function(x) grep(x, names(RW[RWi])))))==0) stop(paste0("I cannot find covariate '", names(RW[RWi]), "' in the survival dataset, \n while it is needed in the requested model (if '", names(RW[RWi]), "' is a \n time-dependent covariate, it needs to be included in the \n survival submodel manually to be properly handled through a \n decomposition of the followup with right censoring and left truncation)."))
            YS_data[[RW[RWi]+1]] <- dataSurv[, which(sapply(sapply(colnames(dataSurv), function(x) grep(x, names(RW[RWi]))), length)!=0)]
          }
        }
      }
      if(TRUE %in% (YS_assoc %in% c("SRE","SRE_ind", "CV", "CS", "CS2", "CV_CS"))){
        # add covariates from the random effects part shared and not already included
        RE <- findbars(formLong[[k]])
        RE_split <- gsub("\\s", "", strsplit(as.character(RE), split=c("\\|"))[[1]])
        RE_elements <- gsub("\\s", "", strsplit(RE_split[[1]], split=c("\\+"))[[1]])
        if(length(which(RE_elements==1))>0) RE_elements[which(RE_elements==1)] <- "Intercept"
        RE_form <- formula(paste(RE_split[2], "~", "-1+", RE_split[1]))
        RE_mat <- as.data.frame(model.matrix(RE_form, model.frame(RE_form, LSurvdat[which(LSurvdat[[id]] %in% dataSurv[[id]]),], na.action=na.pass)))
        colnames(RE_mat) <- RE_elements
        removeVar <- NULL
        for(rmtvar in 1:length(strsplit(colnames(RE_mat), ":"))){ # remove any component that contains a timeVar because it will not be useful
          if(TRUE %in% (c(1, timeVar, c(paste0("f", 1:NFT, "(", timeVar, ")"))) %in% strsplit(colnames(RE_mat), ":")[[rmtvar]]) |
             TRUE %in% (names(YS_data) %in% strsplit(colnames(RE_mat), ":")[[rmtvar]])) removeVar <- c(removeVar, rmtvar)
        }
        colNvar <- colnames(RE_mat)
        RE_mat <- RE_mat[,-c(which(colNvar %in% c("Intercept")), removeVar)]
        if(is.null(dim(RE_mat))){
          if(!(colNvar[-c(which(colNvar %in% c("Intercept")), removeVar)] %in% names(YS_data))){
            YS_data <- append(YS_data, list(RE_mat))
            names(YS_data)[length(names(YS_data))] <- colNvar[-c(which(colNvar %in% c("Intercept")), removeVar)]
          }
        }else if(dim(RE_mat)[2]>0){
          names(RE_mat) <- gsub(":", ".X.", gsub("\\s", ".", names(RE_mat)))
          YS_data <- append(YS_data, RE_mat)
        }
      }
      if(YS_assoc[k]%in%c("CV", "CS", "CS2")){ # one vector
        assign(paste0(YS_assoc[k], "_L", k, "_S", m), c(as.integer(dataSurv[,id]))) # unique id set up after cox expansion
        YS_data <- append(YS_data, list(get(paste0(YS_assoc[k], "_L", k, "_S", m))))
        names(YS_data)[length(names(YS_data))] <- paste0(YS_assoc[k], "_L", k, "_S", m)
      }else if(YS_assoc[k]%in%c("SRE")){
        assign(paste0(YS_assoc[k], "_L", k, "_S", m), c(as.integer(dataSurv[,id]))) # unique id set up after cox expansion
        YS_data <- append(YS_data, list(get(paste0(YS_assoc[k], "_L", k, "_S", m))))
        names(YS_data)[length(names(YS_data))] <- paste0(YS_assoc[k], "_L", k, "_S", m)
      }else if(YS_assoc[k]%in%c("SRE_ind")){
        RE <- findbars(formLong[[k]])
        RE_split <- gsub("\\s", "", strsplit(as.character(RE), split=c("\\|"))[[1]])
        RE_elements <- gsub("\\s", "", strsplit(RE_split[[1]], split=c("\\+"))[[1]])
        if(length(which(RE_elements==1))>0) RE_elements[which(RE_elements==1)] <- "Intercept"
        for(i in 1:length(RE_elements)){
          assign(paste0("SRE_", RE_elements[i], "_L", k, "_S", m), length(unique(c(as.integer(dataSurv[,id]))))*(i-1+CLid) + c(as.integer(dataSurv[,id]))) # unique id set up after cox expansion
          YS_data <- append(YS_data, list(get(paste0("SRE_", RE_elements[i], "_L", k, "_S", m))))
          names(YS_data)[length(names(YS_data))] <- paste0("SRE_", RE_elements[i], "_L", k, "_S", m)
        }
        if(corLong) CLid <- CLid + length(RE_elements)
      }else if(YS_assoc[k]%in%c("CV_CS")){ # need two vectors for current value and current slope
        assign(paste0("CV_L", k, "_S", m), c(as.integer(dataSurv[,id])))
        YS_data <- append(YS_data, list(get(paste0("CV_L", k, "_S", m))))
        names(YS_data)[length(names(YS_data))] <- paste0("CV_L", k, "_S", m)
        assign(paste0("CS_L", k, "_S", m), c(as.integer(dataSurv[,id])))
        YS_data <- append(YS_data, list(get(paste0("CS_L", k, "_S", m))))
        names(YS_data)[length(names(YS_data))] <- paste0("CS_L", k, "_S", m)
      }else{ # SRE_ind => copy (not done yet)

      }
    }
  }
  if(!is.null(RES)){
    RE_splitS <- gsub("\\s", "", strsplit(as.character(RES), split=c("\\|"))[[1]])
    RE_elementS <- gsub("\\s", "", strsplit(RE_splitS[[1]], split=c("\\+"))[[1]])
    RE_formS <- formula(paste(RE_splitS[2], "~", "-1+", RE_splitS[1]))
    RE_matS <- model.matrix(RE_formS, model.frame(RE_formS, dataSurv, na.action=na.pass))
    if(length(which(RE_elementS==1))>0) RE_elementS[which(RE_elementS==1)] <- "Intercept"
    colnames(RE_matS) <- RE_elementS
    colnames(RE_matS) <- gsub("\\s", ".", colnames(RE_matS))
    colnames(RE_matS) <- sub("\\(","", colnames(RE_matS))
    colnames(RE_matS) <- sub(")","", colnames(RE_matS))
  }else RE_matS <- NULL
  if(!is.null(id)){
    if(!id %in% names(YS_data)){ # include id for random effect if not already included for association
      YS_data <- append(YS_data, list(dataSurv[,id]))
      names(YS_data)[length(names(YS_data))] <- id
    }
  }
  names(YS_data) <- sub("\\(","", names(YS_data))
  names(YS_data) <- sub(")","", names(YS_data))
  return(list(YS_data=YS_data, YSformF=YSformF, RE_matS=RE_matS))
}

#' Setup outcome for longitudinal marker
#' @description Setup outcome for longitudinal marker  (internal function)
#' input:
#' @param formula with lme4 format (fixed effects and random effects in the same object)
#' @param dataset that contains the outcome
#' @param family of the outcome (given to check if the distribution matches but the check is not done yet)
#' @param k identifies the longitudinal marker among 1:K markers
#' output:
#' @return YL.name name of the outcome
#' @return YL values of the outcome
setup_Y_model <- function(formula, dataset, family, k){
  dataF <- model.frame(subbars(formula), dataset, na.action=NULL)
  YL.name <- paste0(as.character(subbars(formula))[2], "_L", k)
  YL <- as.vector(model.response(dataF))
  # check if y distribution matches with family here?
  return(list(YL.name, YL))
}

#' Setup fixed effects part for longitudinal marker k
#' @description Setup fixed effects part for longitudinal marker k  (internal function)
#' input:
#' @param formula with lme4 format (fixed effects and random effects in the same object)
#' @param dataset that contains the outcome
#' @param timeVar name of time variable
#' @param k identifies the longitudinal marker among 1:K markers
#'@param dataOnly boolean for internal use, indicates if only preparing data (i.e., not fitting the model)
#' output:
#' @return colnames(FE) names of the fixed effects (interactions are separated
#' by ".X." instead of ":" to facilitate their manipulation)
#' @return FE values of the fixed effects
setup_FE_model <- function(formula, dataset, timeVar, k, dataOnly){
  FE_form <- nobars(formula)
  # if(length(which(sapply(dataset, class)=="factor"))>0 & !dataOnly){ # deal with factors when modalities are missing
  #   factors_columns <- which(sapply(dataset, class)=="factor")
  #   for(fctc in factors_columns){
  #     if(length(unique(dataset[, fctc]))< length(levels(dataset[, fctc]))) dataset[, fctc] <- as.integer(dataset[, fctc])-1
  #   }
  # }
  FE <- model.matrix(FE_form, model.frame(FE_form, dataset, na.action=na.pass))

  # Fix factor handling when no intercept is present (-1 in formula)
  # Check if formula has no intercept
  has_intercept <- attr(terms(FE_form), "intercept") == 1

  if (!has_intercept) {
    # Get original data frame to identify factor variables
    mf <- model.frame(FE_form, dataset, na.action=na.pass)

    # Find factor variables in the original data
    factor_vars <- names(mf)[sapply(mf, is.factor)]

    # For each factor variable, check if all levels are included in model matrix
    for (fvar in factor_vars) {
      if (fvar %in% all.vars(FE_form)) {
        factor_levels <- levels(mf[[fvar]])

        # Only handle binary factors (most common case)
        if (length(factor_levels) == 2) {
          # Check if both factor levels appear as standalone columns (not in interactions)
          factor_col_0 <- paste0(fvar, factor_levels[1])
          factor_col_1 <- paste0(fvar, factor_levels[2])

          if (factor_col_0 %in% colnames(FE) && factor_col_1 %in% colnames(FE)) {
            # Check that these are main effect columns, not part of interactions
            # Main effect columns should have simple names like "trt0", "trt1"
            main_effect_cols <- c(factor_col_0, factor_col_1)

            # Find which columns are main effects vs interactions
            all_cols <- colnames(FE)
            main_cols <- all_cols[!grepl(":", all_cols)]

            # Only modify if both levels appear as main effects
            if (all(main_effect_cols %in% main_cols)) {
              # Remove the first level column (reference level)
              cols_to_keep <- colnames(FE) != factor_col_0
              FE <- FE[, cols_to_keep, drop = FALSE]
            }
          }
        }
      }
    }
  }

  #if(colnames(FE)[1]=="(Intercept)") colnames(FE)[1] <- "Intercept"
  colnames(FE) <- gsub(":", ".X.", gsub("\\s", ".", colnames(FE)))
  colnames(FE) <- sub("\\(","", colnames(FE))
  colnames(FE) <- sub(")","", colnames(FE))
  OFS <- model.offset(model.frame(FE_form, dataset, na.action=na.pass))
  if(is.null(OFS)){
    OFS <- rep(0, dim(FE)[1])
  }
  return(list(colnames(FE), FE, OFS))
}

#' Setup random effects part for longitudinal marker k
#' @description Setup random effects part for longitudinal marker k  (internal function)
#' input:
#' @param formula with lme4 format (fixed effects and random effects in the same object)
#' @param dataset that contains the outcome
#' @param k identifies the longitudinal marker among 1:K markers
#' output:
#' @return colnames(RE_mat) names of the random effects
#' @return RE_mat values of the random effects
setup_RE_model <- function(formula, dataset, k){
  RE <- findbars(formula)
  if(length(RE)==0) stop("No random effects found, the longitudinal part must at least contain one random effect per marker.")
  RE_split <- gsub("\\s", "", strsplit(as.character(RE), split=c("\\|"))[[1]])
  RE_elements <- gsub("\\s", "", strsplit(RE_split[[1]], split=c("\\+"))[[1]])
  RE_form <- formula(paste(RE_split[2], "~", "-1+", RE_split[1]))
  RE_mat <- model.matrix(RE_form, model.frame(RE_form, dataset, na.action=na.pass))
  if(length(which(RE_elements==1))>0) RE_elements[which(RE_elements==1)] <- "Intercept"
  colnames(RE_mat) <- RE_elements
  #if(colnames(RE_mat)[1]=="(Intercept)") colnames(RE_mat)[1] <- "Intercept"
  colnames(RE_mat) <- gsub("\\s", ".", colnames(RE_mat))
  colnames(RE_mat) <- sub("\\(","", colnames(RE_mat))
  colnames(RE_mat) <- sub(")","", colnames(RE_mat))
  return(list(colnames(RE_mat), RE_mat))
}

#' Setup scopy
#' @description Setup weights for non-linear effects
#' input:
#' @param n number of knots
#' output:
#' @return Matrix with weights for scopy
INLAjoint.scopy.define <- function(n = 5L)
{
  stopifnot(n == 2 || n >= 5)
  if (n == 2) {
    VV = cbind(1, c(-0.5, 0.5))
  }
  else {
    Q <- INLAjoint.rw2(n, scale.model = TRUE)
    e <- eigen(Q)
    idx.remove <- (n - 1):n
    V <- e$vectors[, -idx.remove, drop = FALSE]
    lambda <- e$values[-idx.remove]
    VV <- matrix(0, n, n)
    VV[, 1:2] <- cbind(rep(1, n), seq(-0.5, 0.5, len = n))
    m <- n - 2
    for (i in seq_len(m)) {
      j <- m - i + 3
      VV[, j] <- V[, i] * sqrt(1/lambda[i])
    }
  }
  return(list(n = n, W = VV))
}

#' Setup rw2
#' @description Setup random walk of order 2
#' input:
#' @param n number of knots
#' @param order order 2
#' @param ... additional arguments passed to INLAjoint.rw
#' output:
#' @return random walk 2
INLAjoint.rw2 <- function (n, order=2L, ...)
{
  INLAjoint.rw(n, order = order, ...)
}

#' Setup rw
#' @description Setup random walk
#' input:
#' @param n number of knots
#' @param order order 1
#' @param sparse boolean, sparsity
#' @param scale.model boolean, scale
#' @param cyclic boolean, cyclic
#' output:
#' @return random walk 1
INLAjoint.rw <- function (n, order = 1L, sparse = TRUE, scale.model = FALSE,
          cyclic = FALSE)
{
  stopifnot(n >= 1L + 2L * order)
  if (scale.model) {
    if (!cyclic) {
      rd <- order
    }
    else {
      if (order == 1L) {
        rd <- 1L
      }
      else {
        rd <- order - 1L
      }
    }
    Q <- INLAjoint.rw(n, order = order, sparse = sparse, scale.model = FALSE,
                 cyclic = cyclic)
    fac <- exp(mean(log(diag(INLAjoint.ginv(as.matrix(Q), rankdef = rd)))))
    Q <- fac * Q
  }
  else {
    if (!cyclic) {
      U <- diff(diag(n), diff = order)
      Q <- t(U) %*% U
    }
    else {
      m <- 2L * order + 1L
      k <- 1L + order
      U <- diff(diag(m), diff = order)
      U <- t(U) %*% U
      Q <- toeplitz(c(U[k, k:m], rep(0, n - m), U[k, m:(k +
                                                          1L)]))
    }
  }
  return(if (sparse) INLA::inla.as.sparse(Q) else Q)
}

#' Setup ginv
#' @description Setup random walk
#' input:
#' @param x input
#' @param tol tolerance
#' @param rankdef rank
#' output:
#' @return ginv
INLAjoint.ginv <- function (x, tol = sqrt(.Machine$double.eps), rankdef = NULL)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x))) {
    stop("'x' must be a numeric or complex matrix")
  }
  xsvd <- svd(x)
  if (is.complex(x)) {
    xsvd$u <- Conj(xsvd$u)
  }
  if (is.null(rankdef) || rankdef == 0L) {
    Positive <- xsvd$d > max(tol * xsvd$d[1L], 0)
  }
  else {
    n <- length(xsvd$d)
    stopifnot(rankdef >= 1L && rankdef <= n)
    Positive <- c(rep(TRUE, n - rankdef), rep(FALSE, rankdef))
  }
  if (all(Positive)) {
    xsvd$v %*% (1/xsvd$d * t(xsvd$u))
  }
  else if (!any(Positive)) {
    array(0, dim(x)[2L:1L])
  }
  else {
    xsvd$v[, Positive, drop = FALSE] %*% ((1/xsvd$d[Positive]) *
                                            t(xsvd$u[, Positive, drop = FALSE]))
  }
}



# RW2
parse_rw2 <- function(formula) {
  form_str <- deparse(formula, width.cutoff = 500L)
  form_str <- paste(form_str, collapse = " ")

  rw2_pattern <- "RW2\\s*\\(([^,]+),\\s*group\\s*=\\s*([^)]+)\\)"
  has_rw2 <- grepl(rw2_pattern, form_str, perl = TRUE)

  if (!has_rw2) {
    return(list(has_rw2 = FALSE, time_var = NULL, group_expr = NULL,
                clean_formula = formula))
  }

  matches <- regmatches(form_str, regexec(rw2_pattern, form_str, perl = TRUE))[[1]]
  time_var <- trimws(matches[2])
  group_expr <- trimws(matches[3])

  # Replace RW2(TIME, group=...) with just TIME
  # This keeps TIME in the formula for initial data processing
  clean_str <- gsub(rw2_pattern, time_var, form_str, perl = TRUE)
  clean_str <- gsub("\\s+", " ", clean_str)
  clean_str <- trimws(clean_str)

  clean_formula <- as.formula(clean_str)

  return(list(
    has_rw2 = TRUE,
    time_var = time_var,
    group_expr = group_expr,
    clean_formula = clean_formula
  ))
}



